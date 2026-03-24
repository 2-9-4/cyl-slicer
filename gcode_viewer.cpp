// gcode_viewer.cpp — Volumetric OpenGL G-code viewer for cylindrical printer
// Lines rendered as 3D beads (layer_height × nozzle_size rectangular cross-section)
// Controls: left-drag=orbit  right/mid-drag=pan  scroll=zoom  R=reset  Q/Esc=quit

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <algorithm>
#include <cmath>
#include <climits>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// ── Slicer parameters (parsed from G-code header; these are fallback defaults) ─

static const float kDefaultLayerHeight = 0.2f;
static const float kDefaultLineWidth   = 0.4f;

// ── Section types ─────────────────────────────────────────────────────────────

enum SType {
    ST_EXTERNAL = 0, ST_INTERNAL, ST_SOLID_INFILL,
    ST_RADIAL_FULL,  ST_RADIAL,   ST_AXIAL,
    ST_SUPPORT,      ST_TRAVEL,   ST_RETRACT, ST_UNKNOWN, ST_COUNT
};

static const char* kNames[ST_COUNT] = {
    "External Perimeter", "Internal Perimeter", "Solid Infill",
    "Infill Radial Full", "Infill Radial",      "Infill Axial",
    "Support",            "Travel",             "Retraction",  "Unknown",
};

// RGB in [0,1]
static const float kColors[ST_COUNT][3] = {
    {0.000f, 0.667f, 0.000f},  // EXTERNAL      #00AA00
    {0.000f, 0.400f, 0.800f},  // INTERNAL      #0066CC
    {0.800f, 0.400f, 0.000f},  // SOLID_INFILL  #CC6600
    {0.800f, 0.000f, 0.800f},  // RADIAL_FULL   #CC00CC
    {1.000f, 0.400f, 0.000f},  // RADIAL        #FF6600
    {0.400f, 0.000f, 0.800f},  // AXIAL         #6600CC
    {0.000f, 0.800f, 0.800f},  // SUPPORT       #00CCCC cyan
    {1.000f, 0.200f, 0.200f},  // TRAVEL        red
    {1.000f, 0.900f, 0.000f},  // RETRACT       yellow
    {0.533f, 0.533f, 0.533f},  // UNKNOWN       #888888
};

static const float kAlpha[ST_COUNT] = {
    1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.10f, 1.0f, 0.6f
};

SType parse_stype(const std::string& s) {
    if (s == "EXTERNAL_PERIMETER") return ST_EXTERNAL;
    if (s == "INTERNAL_PERIMETER") return ST_INTERNAL;
    if (s == "SOLID_INFILL")       return ST_SOLID_INFILL;
    if (s == "INFILL_RADIAL_FULL") return ST_RADIAL_FULL;
    if (s == "INFILL_RADIAL")      return ST_RADIAL;
    if (s == "INFILL_AXIAL")       return ST_AXIAL;
    if (s == "SUPPORT")            return ST_SUPPORT;
    return ST_UNKNOWN;
}

// ── Vec3 helpers ──────────────────────────────────────────────────────────────

struct Vec3 { float x, y, z; };

inline Vec3 v3sub(Vec3 a, Vec3 b)   { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline Vec3 v3add(Vec3 a, Vec3 b)   { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline Vec3 v3scale(Vec3 a, float s){ return {a.x*s,   a.y*s,   a.z*s};   }
inline Vec3 v3neg(Vec3 a)           { return {-a.x,    -a.y,    -a.z};     }
inline float v3dot(Vec3 a, Vec3 b)  { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline Vec3 v3cross(Vec3 a, Vec3 b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
inline Vec3 v3norm(Vec3 a) {
    float l = sqrtf(v3dot(a, a));
    return l > 1e-8f ? v3scale(a, 1.f / l) : Vec3{0, 1, 0};
}

Vec3 cyl_to_cart(float theta_deg, float ax, float r) {
    float t = std::fmod(theta_deg, 360.f) * (float)M_PI / 180.f;
    return { ax, r * std::cos(t), r * std::sin(t) };
}

// Radial direction outward from cylinder axis (X-axis)
static Vec3 radial_dir(Vec3 p) {
    float len = sqrtf(p.y * p.y + p.z * p.z);
    return len > 1e-6f ? Vec3{0, p.y / len, p.z / len} : Vec3{0, 1, 0};
}

// ── Vertex with position + normal ─────────────────────────────────────────────

struct Vert { float px, py, pz, nx, ny, nz; };

// ── Tube mesh generation ──────────────────────────────────────────────────────
//
// Sweeps a rectangle (lh × lw) along the centerline path.
// lh = layer_height → radial extent (in/out of cylinder surface)
// lw = line_width   → transverse extent (along cylinder surface, ⊥ path)
//
// Cross-section corner convention:
//   c[0] = outer + side  (+r, +s)
//   c[1] = outer - side  (+r, -s)
//   c[2] = inner - side  (-r, -s)
//   c[3] = inner + side  (-r, +s)

static void build_tube(
    std::vector<Vert>&     verts,
    std::vector<uint32_t>& idx,
    const std::vector<Vec3>& path,
    float lh, float lw)
{
    int N = (int)path.size();
    if (N < 2) return;

    float hw = lw * 0.5f;
    float hh = lh * 0.5f;

    // Tangent at each vertex
    std::vector<Vec3> tang(N);
    for (int i = 0; i < N; ++i) {
        if (i == 0)
            tang[i] = v3norm(v3sub(path[1], path[0]));
        else if (i == N - 1)
            tang[i] = v3norm(v3sub(path[N-1], path[N-2]));
        else
            tang[i] = v3norm(v3add(v3norm(v3sub(path[i],   path[i-1])),
                                   v3norm(v3sub(path[i+1], path[i]))));
    }

    // Cross-section (4 corners + r/s axes) at each vertex
    struct XS { Vec3 c[4]; Vec3 r, s; };
    std::vector<XS> xs(N);

    // s = cross(r, tang): always perpendicular to both radial and path direction.
    // For circumferential paths: cross(radial, circumferential) = axial  ✓
    // For axial paths:           cross(radial, axial)           = circumferential ✓
    // No path-type classification needed.
    for (int i = 0; i < N; ++i) {
        Vec3 r = radial_dir(path[i]);
        Vec3 s = v3norm(v3cross(r, tang[i]));
        if (v3dot(s, s) < 0.01f) {
            // Degenerate: tang parallel to r — fall back to axial direction
            s = v3norm(v3cross(r, Vec3{1, 0, 0}));
        }
        xs[i].r = r; xs[i].s = s;
    }

    // Compute corners for each cross-section
    for (int i = 0; i < N; ++i) {
        Vec3 p = path[i];
        xs[i].c[0] = {p.x + hh*xs[i].r.x + hw*xs[i].s.x, p.y + hh*xs[i].r.y + hw*xs[i].s.y, p.z + hh*xs[i].r.z + hw*xs[i].s.z};
        xs[i].c[1] = {p.x + hh*xs[i].r.x - hw*xs[i].s.x, p.y + hh*xs[i].r.y - hw*xs[i].s.y, p.z + hh*xs[i].r.z - hw*xs[i].s.z};
        xs[i].c[2] = {p.x - hh*xs[i].r.x - hw*xs[i].s.x, p.y - hh*xs[i].r.y - hw*xs[i].s.y, p.z - hh*xs[i].r.z - hw*xs[i].s.z};
        xs[i].c[3] = {p.x - hh*xs[i].r.x + hw*xs[i].s.x, p.y - hh*xs[i].r.y + hw*xs[i].s.y, p.z - hh*xs[i].r.z + hw*xs[i].s.z};
    }

    // Emit a quad (a,b,c,d in CCW order) with the given face normal
    auto add_quad = [&](Vec3 a, Vec3 b, Vec3 c, Vec3 d, Vec3 n) {
        uint32_t base = (uint32_t)verts.size();
        verts.push_back({a.x, a.y, a.z, n.x, n.y, n.z});
        verts.push_back({b.x, b.y, b.z, n.x, n.y, n.z});
        verts.push_back({c.x, c.y, c.z, n.x, n.y, n.z});
        verts.push_back({d.x, d.y, d.z, n.x, n.y, n.z});
        idx.push_back(base);     idx.push_back(base + 1); idx.push_back(base + 2);
        idx.push_back(base);     idx.push_back(base + 2); idx.push_back(base + 3);
    };

    // Lateral faces between consecutive cross-sections
    for (int i = 0; i < N - 1; ++i) {
        const XS& x0 = xs[i], &x1 = xs[i + 1];
        Vec3 r = v3norm(v3add(x0.r, x1.r));
        Vec3 s = v3norm(v3add(x0.s, x1.s));

        // Outer face (+r normal)
        add_quad(x0.c[0], x0.c[1], x1.c[1], x1.c[0], r);
        // Inner face (-r normal)
        add_quad(x0.c[2], x0.c[3], x1.c[3], x1.c[2], v3neg(r));
        // +side face (+s normal)
        add_quad(x0.c[3], x0.c[0], x1.c[0], x1.c[3], s);
        // -side face (-s normal)
        add_quad(x0.c[1], x0.c[2], x1.c[2], x1.c[1], v3neg(s));
    }

    // Start cap (normal = -tang[0])
    {
        const XS& x = xs[0];
        Vec3 n = v3neg(tang[0]);
        add_quad(x.c[0], x.c[3], x.c[2], x.c[1], n);
    }
    // End cap (normal = +tang[N-1])
    {
        const XS& x = xs[N - 1];
        add_quad(x.c[0], x.c[1], x.c[2], x.c[3], tang[N - 1]);
    }
}

// ── Octahedron marker (retraction point) ─────────────────────────────────────

static void build_sphere_marker(
    std::vector<Vert>&     verts,
    std::vector<uint32_t>& idx,
    Vec3 center, float r)
{
    Vec3 pts[6] = {
        {center.x + r, center.y,     center.z    },
        {center.x - r, center.y,     center.z    },
        {center.x,     center.y + r, center.z    },
        {center.x,     center.y - r, center.z    },
        {center.x,     center.y,     center.z + r},
        {center.x,     center.y,     center.z - r},
    };
    int faces[8][3] = {
        {0,2,4}, {0,4,3}, {0,3,5}, {0,5,2},
        {1,4,2}, {1,3,4}, {1,5,3}, {1,2,5},
    };
    for (auto& f : faces) {
        Vec3 a = pts[f[0]], b = pts[f[1]], c = pts[f[2]];
        Vec3 n = v3norm(v3cross(v3sub(b, a), v3sub(c, a)));
        uint32_t base = (uint32_t)verts.size();
        verts.push_back({a.x, a.y, a.z, n.x, n.y, n.z});
        verts.push_back({b.x, b.y, b.z, n.x, n.y, n.z});
        verts.push_back({c.x, c.y, c.z, n.x, n.y, n.z});
        idx.push_back(base); idx.push_back(base + 1); idx.push_back(base + 2);
    }
}

// ── Cone/arrow marker (playback position indicator) ──────────────────────────

static void build_cone_marker(
    std::vector<Vert>&     verts,
    std::vector<uint32_t>& idx,
    Vec3 tip, Vec3 direction, float length, float radius)
{
    // Cone pointing in direction, with tip at 'tip'
    Vec3 dir = v3norm(direction);
    Vec3 base_center = {tip.x - dir.x * length, tip.y - dir.y * length, tip.z - dir.z * length};
    
    // Build orthonormal basis
    Vec3 up = fabsf(dir.x) < 0.9f ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
    Vec3 right = v3norm(v3cross(dir, up));
    Vec3 fwd = v3norm(v3cross(right, dir));
    
    int slices = 16;
    uint32_t base_idx = (uint32_t)verts.size();
    
    // Tip vertex
    verts.push_back({tip.x, tip.y, tip.z, dir.x, dir.y, dir.z});
    
    // Base circle vertices
    for (int i = 0; i <= slices; ++i) {
        float angle = 2.f * (float)M_PI * i / slices;
        float c = cosf(angle), s = sinf(angle);
        Vec3 p = {
            base_center.x + radius * (c * right.x + s * fwd.x),
            base_center.y + radius * (c * right.y + s * fwd.y),
            base_center.z + radius * (c * right.z + s * fwd.z)
        };
        Vec3 n = v3norm(v3sub(p, base_center));
        verts.push_back({p.x, p.y, p.z, n.x, n.y, n.z});
    }
    
    // Cone sides (triangles from tip to base edge)
    for (int i = 0; i < slices; ++i) {
        idx.push_back(base_idx);
        idx.push_back(base_idx + 1 + i);
        idx.push_back(base_idx + 2 + i);
    }
    
    // Base cap (fan from base center)
    uint32_t cap_base = (uint32_t)verts.size();
    Vec3 cap_n = v3neg(dir);
    verts.push_back({base_center.x, base_center.y, base_center.z, cap_n.x, cap_n.y, cap_n.z});
    for (int i = 0; i <= slices; ++i) {
        float angle = 2.f * (float)M_PI * i / slices;
        float c = cosf(angle), s = sinf(angle);
        Vec3 p = {
            base_center.x + radius * (c * right.x + s * fwd.x),
            base_center.y + radius * (c * right.y + s * fwd.y),
            base_center.z + radius * (c * right.z + s * fwd.z)
        };
        verts.push_back({p.x, p.y, p.z, cap_n.x, cap_n.y, cap_n.z});
    }
    for (int i = 0; i < slices; ++i) {
        idx.push_back(cap_base);
        idx.push_back(cap_base + 2 + i);
        idx.push_back(cap_base + 1 + i);
    }
}

// ── Path group ────────────────────────────────────────────────────────────────

struct LayerRange {
    int      layer_num;
    int      seq;        // sequence within layer (-1 if not applicable)
    uint32_t idx_offset;
    uint32_t idx_count;
};

struct PathGroup {
    SType type;
    std::vector<Vert>     verts;
    std::vector<uint32_t> idx;
    std::vector<LayerRange> layer_ranges;
    GLuint vao = 0, vbo = 0, ebo = 0;

    // Build state
    std::vector<Vec3> cur_path;
    float tube_lh = kDefaultLayerHeight;
    float tube_lw = kDefaultLineWidth;

    void add_pt(Vec3 p)  { cur_path.push_back(p); }

    std::function<void(SType, uint32_t, uint32_t, Vec3)> on_flush;

    void flush_path() {
        uint32_t before = (uint32_t)idx.size();
        Vec3 end_pos = cur_path.size() >= 1 ? cur_path.back() : Vec3{0,0,0};
        if (cur_path.size() >= 2)
            build_tube(verts, idx, cur_path, tube_lh, tube_lw);
        cur_path.clear();
        uint32_t after = (uint32_t)idx.size();
        if (after > before && on_flush)
            on_flush(type, before, after - before, end_pos);
    }

    void break_path() { flush_path(); }
};

struct DrawCmd {
    SType    type;
    int      layer;
    int      seq;        // sequence within layer (0-based)
    uint32_t idx_offset;
    uint32_t idx_count;
    Vec3     end_pos;    // end position of this path (for playback pointer)
    std::vector<std::string> lines;  // raw G-code lines that make up this path
};

struct GcodeData {
    PathGroup groups[ST_COUNT];
    struct RetractPt { Vec3 pos; int layer; int seq; };
    std::vector<RetractPt>        retract_pts;  // (position, layer_num, seq)
    std::vector<DrawCmd>              draw_cmds;
    std::unordered_map<int, int>      layer_seq_count; // layer -> # seqs
    float x_min = 1e9f, x_max = -1e9f;
    float r_min = 1e9f, r_max = -1e9f;
    int   layer_min = INT_MAX, layer_max = 0;
};

// ── G-code parsing ────────────────────────────────────────────────────────────

GcodeData load_gcode(const char* path) {
    GcodeData d;
    for (int i = 0; i < ST_COUNT; ++i) d.groups[i].type = (SType)i;

    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open: " << path << "\n"; return d; }

    float cur_x = 0.f, cur_t = 0.f, cur_e = 0.f;
    float radius = 10.f;
    bool  has_pos = false;
    bool  in_ext  = false;
    SType ext_type = ST_UNKNOWN;
    SType sec_type = ST_UNKNOWN;
    int   cur_layer = 0;
    int   cur_layer_seq = 0;

    // Parameters parsed from G-code header comments
    float lh = kDefaultLayerHeight, lw = kDefaultLineWidth;

    uint32_t layer_idx_start[ST_COUNT] = {};

    // Per-type line buffers: lines accumulated since last flush for each section type
    std::vector<std::string> pending_lines[ST_COUNT];

    // Record each flushed path as a DrawCmd ordered by sequence within its layer
    for (int i = 0; i < ST_COUNT; ++i) {
        d.groups[i].on_flush = [&, i](SType t, uint32_t off, uint32_t cnt, Vec3 end_pos) {
            DrawCmd cmd{t, cur_layer, cur_layer_seq++, off, cnt, end_pos, {}};
            cmd.lines = std::move(pending_lines[i]);
            pending_lines[i].clear();
            d.draw_cmds.push_back(std::move(cmd));
        };
    }

    auto close_layer = [&]() {
        // Flush in-progress paths before snapshotting layer ranges
        for (int i = 0; i < ST_COUNT; ++i)
            d.groups[i].flush_path();
        for (int i = 0; i < ST_COUNT; ++i) {
            uint32_t start = layer_idx_start[i];
            uint32_t end   = (uint32_t)d.groups[i].idx.size();
            if (end > start)
                d.groups[i].layer_ranges.push_back({cur_layer, -1, start, end - start});
            layer_idx_start[i] = end;
        }
        cur_layer_seq = 0;  // reset for next layer
    };

    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        if (line[0] == ';') {
            if (line.rfind(";LAYER:", 0) == 0) {
                close_layer();
                std::string layer_str = line.substr(7);
                if (!layer_str.empty() && (std::isdigit(layer_str[0]) || layer_str[0] == '-')) {
                    cur_layer = std::stoi(layer_str);
                    d.layer_min = std::min(d.layer_min, cur_layer);
                    d.layer_max = std::max(d.layer_max, cur_layer);
                }
                in_ext = false;
            } else if (line.rfind(";LAYER_HEIGHT:", 0) == 0) {
                lh = std::stof(line.substr(14));
                for (auto& grp : d.groups) grp.tube_lh = lh;
            } else if (line.rfind(";LINE_WIDTH:", 0) == 0) {
                lw = std::stof(line.substr(12));
                for (auto& grp : d.groups) grp.tube_lw = lw;
            } else if (line.rfind(";Z:", 0) == 0) {
                radius = std::stof(line.substr(3));
                d.r_min = std::min(d.r_min, radius);
                d.r_max = std::max(d.r_max, radius);
            } else if (line.rfind(";TYPE:", 0) == 0) {
                SType nt = parse_stype(line.substr(6));
                if (nt != sec_type) { sec_type = nt; in_ext = false; }
            }
            continue;
        }

        std::istringstream ss(line);
        std::string cmd; ss >> cmd;
        if (cmd != "G0" && cmd != "G1") continue;

        bool rapid = (cmd == "G0");
        float nx = cur_x, nt = cur_t, ne = cur_e;
        bool has_x_tok = false, has_t_tok = false, has_e_tok = false;

        std::string tok;
        while (ss >> tok) {
            if (tok.size() < 2) continue;
            try {
                float v = std::stof(tok.substr(1));
                if      (tok[0] == 'X') { nx = v; has_x_tok = true; }
                else if (tok[0] == 'Y') { nt = v; has_t_tok = true; }
                else if (tok[0] == 'E') { ne = v; has_e_tok = true; }
            } catch (...) {}
        }

        // Detect retraction: G1 with E going down and no X/Y motion
        if (!rapid && has_e_tok && ne < cur_e - 1e-6f && !has_x_tok && !has_t_tok && has_pos)
            d.retract_pts.push_back({cyl_to_cart(cur_t, cur_x, radius), cur_layer, cur_layer_seq});
        cur_e = ne;

        d.x_min = std::min(d.x_min, nx);
        d.x_max = std::max(d.x_max, nx);
        Vec3 nv = cyl_to_cart(nt, nx, radius);

        if (rapid) {
            in_ext = false;
            // End any ongoing extrusion path
            for (int i = 0; i < ST_COUNT; ++i)
                if (i != ST_TRAVEL) d.groups[i].flush_path();

            if (has_pos) {
                auto& g = d.groups[ST_TRAVEL];
                g.flush_path();  // flush any previous (should be empty after immediate-flush below)
                pending_lines[ST_TRAVEL].push_back(line);
                g.add_pt(cyl_to_cart(cur_t, cur_x, radius));
                g.add_pt(nv);
                g.flush_path();  // flush immediately so seq is assigned before the next G1
            }
        } else {
            // Compute angular delta with wrap-around handling
            float raw_delta_t = nt - cur_t;
            // For tessellation, we want the actual angular travel, not shortest path
            // A move from Y360 to Y720 is a full 360° circle, not 0°
            // But we still need to handle cases where the slicer uses Y0 vs Y360 interchangeably
            
            // Check if this is a "large" angular move (> 180°) - likely intentional full/partial circle
            bool large_angular = fabsf(raw_delta_t) > 180.f;
            
            // Compute Cartesian positions
            Vec3 cur_pos = cyl_to_cart(cur_t, cur_x, radius);
            Vec3 new_pos = cyl_to_cart(nt, nx, radius);
            float dx = new_pos.x - cur_pos.x;
            float dy = new_pos.y - cur_pos.y;
            float dz = new_pos.z - cur_pos.z;
            float dist_sq = dx*dx + dy*dy + dz*dz;
            
            // Skip only truly degenerate segments (no axial AND no significant angular motion)
            float axial_dist_sq = (nx - cur_x) * (nx - cur_x);
            bool has_axial_motion = axial_dist_sq > 1e-6f;
            bool has_angular_motion = large_angular || (dist_sq > 1e-6f && fabsf(raw_delta_t) > 0.1f);
            
            if (!has_axial_motion && !has_angular_motion) {
                cur_x = nx; cur_t = nt;
                continue;
            }

            if (!in_ext || sec_type != ext_type) {
                // Type change: end previous path in old group
                if (in_ext) d.groups[ext_type].flush_path();
                in_ext   = true;
                ext_type = sec_type;
                auto& g = d.groups[ext_type];
                // Seed new path with current position
                if (has_pos) g.add_pt(cur_pos);
            }
            pending_lines[ext_type].push_back(line);  // record G1 in this path's feed

            // Tessellate angular spans
            float delta_t = raw_delta_t;
            if (!large_angular) {
                // For small moves, use shortest path across 0/360 boundary
                if (delta_t > 180.f)      delta_t -= 360.f;
                else if (delta_t < -180.f) delta_t += 360.f;
            }
            // For large angular moves (>180°), use the raw delta_t as-is (full circle, etc.)
            
            int steps = std::max(1, (int)(std::abs(delta_t) / 2.f));
            for (int s = 1; s <= steps; ++s) {
                float frac = (float)s / (float)steps;
                float it = cur_t + delta_t * frac;
                // Keep interpolated angle in [0, 360) for consistent Cartesian conversion
                it = std::fmod(it, 360.f);
                if (it < 0.f) it += 360.f;
                float ix = cur_x + (nx - cur_x) * frac;
                d.groups[ext_type].add_pt(cyl_to_cart(it, ix, radius));
            }
        }

        cur_x = nx; cur_t = nt; has_pos = true;
    }

    close_layer();  // flush last layer
    if (d.layer_min == INT_MAX) { d.layer_min = 0; d.layer_max = 0; }
    for (auto& grp : d.groups) grp.on_flush = nullptr;
    return d;
}

// ── Cylinder mesh (occluder) ──────────────────────────────────────────────────

struct Mesh { GLuint vao = 0, vbo = 0, ebo = 0; int count = 0; };

Mesh build_cylinder(float r, float x0, float x1, int slices = 64) {
    // Vertex layout: px py pz nx ny nz  (6 floats)
    std::vector<float>    v;
    std::vector<uint32_t> idx;

    for (int xi = 0; xi < 2; ++xi) {
        float ax = (xi == 0) ? x0 : x1;
        for (int i = 0; i <= slices; ++i) {
            float t = 2.f * (float)M_PI * i / slices;
            float cy = std::cos(t), cz = std::sin(t);
            v.push_back(ax);       v.push_back(r * cy);  v.push_back(r * cz);
            v.push_back(0.f);      v.push_back(cy);      v.push_back(cz);
        }
    }
    for (int i = 0; i <= slices; ++i) {
        idx.push_back((uint32_t)i);
        idx.push_back((uint32_t)(slices + 1 + i));
    }

    Mesh m;
    glGenVertexArrays(1, &m.vao);
    glGenBuffers(1, &m.vbo);
    glGenBuffers(1, &m.ebo);
    glBindVertexArray(m.vao);
    glBindBuffer(GL_ARRAY_BUFFER, m.vbo);
    glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(v.size() * sizeof(float)), v.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, (GLsizeiptr)(idx.size() * sizeof(uint32_t)), idx.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    m.count = (int)idx.size();
    return m;
}

// ── Shaders ───────────────────────────────────────────────────────────────────

static const char* kVert = R"glsl(
#version 330 core
layout(location=0) in vec3 pos;
layout(location=1) in vec3 normal;
uniform mat4 MVP;
out vec3 vNormal;
void main() {
    gl_Position = MVP * vec4(pos, 1.0);
    vNormal = normal;
}
)glsl";

// Simple Lambert shading with a fixed world-space light.
// abs(dot) makes both faces of each quad respond to the light,
// so the viewer is not confused by inside/outside orientation.
static const char* kFrag = R"glsl(
#version 330 core
in  vec3 vNormal;
uniform vec4 color;
uniform vec3 lightDir;   // normalised, set once at startup
out vec4 fragColor;
void main() {
    float diff    = abs(dot(normalize(vNormal), lightDir));
    float ambient = 0.30;
    float lit     = ambient + (1.0 - ambient) * diff;
    fragColor = vec4(color.rgb * lit, color.a);
}
)glsl";

GLuint make_program() {
    auto compile = [](GLenum type, const char* src) {
        GLuint s = glCreateShader(type);
        glShaderSource(s, 1, &src, nullptr);
        glCompileShader(s);
        GLint ok; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
        if (!ok) { char b[512]; glGetShaderInfoLog(s, 512, nullptr, b); std::cerr << b << "\n"; }
        return s;
    };
    GLuint vs = compile(GL_VERTEX_SHADER,   kVert);
    GLuint fs = compile(GL_FRAGMENT_SHADER, kFrag);
    GLuint p  = glCreateProgram();
    glAttachShader(p, vs); glAttachShader(p, fs);
    glLinkProgram(p);
    glDeleteShader(vs); glDeleteShader(fs);
    return p;
}

// ── Camera ────────────────────────────────────────────────────────────────────

struct Camera {
    float az = 0.5f, el = 0.3f, dist = 200.f;
    glm::vec3 target{0, 0, 0};

    glm::vec3 eye() const {
        return target + glm::vec3(
            dist * std::cos(el) * std::sin(az),
            dist * std::sin(el),
            dist * std::cos(el) * std::cos(az));
    }
    glm::mat4 view() const {
        return glm::lookAt(eye(), target, glm::vec3(0, 1, 0));
    }
};

// ── App state ─────────────────────────────────────────────────────────────────

struct App {
    Camera cam, cam0;
    bool   visible[ST_COUNT];
    double mx = 0, my = 0;
    bool   drag = false, pan = false;
    int    W = 1400, H = 900;
    int    show_layer_min = 0, show_layer_max = 0;
    int    total_layer_min = 0, total_layer_max = 0;
    float  show_seq_float = -1.f;  // float for fine control
    int    show_seq = -1;          // -1 = show all; used only when single layer selected
    int    seq_layer = INT_MIN;    // which layer show_seq is tracking (for reset on change)
    int    feed_prev_seq = -2;    // previous show_seq, used to detect changes for auto-scroll
} g;

void cb_mouse(GLFWwindow* w, int btn, int act, int) {
    if (ImGui::GetIO().WantCaptureMouse) return;
    glfwGetCursorPos(w, &g.mx, &g.my);
    if (btn == GLFW_MOUSE_BUTTON_LEFT)                                     g.drag = act == GLFW_PRESS;
    if (btn == GLFW_MOUSE_BUTTON_RIGHT || btn == GLFW_MOUSE_BUTTON_MIDDLE) g.pan  = act == GLFW_PRESS;
}

void cb_pos(GLFWwindow*, double x, double y) {
    double dx = x - g.mx, dy = y - g.my;
    g.mx = x; g.my = y;
    if (ImGui::GetIO().WantCaptureMouse) return;
    if (g.drag) {
        g.cam.az -= (float)(dx * 0.005);
        g.cam.el  = glm::clamp(g.cam.el + (float)(dy * 0.005), -1.5f, 1.5f);
    }
    if (g.pan) {
        glm::vec3 fwd   = glm::normalize(g.cam.target - g.cam.eye());
        glm::vec3 right = glm::normalize(glm::cross(fwd, {0, 1, 0}));
        glm::vec3 up    = glm::cross(right, fwd);
        float s = g.cam.dist * 0.001f;
        g.cam.target -= right * (float)(dx * s);
        g.cam.target += up    * (float)(dy * s);
    }
}

void cb_scroll(GLFWwindow*, double, double dy) {
    if (ImGui::GetIO().WantCaptureMouse) return;
    g.cam.dist = glm::clamp(g.cam.dist * std::pow(0.9f, (float)dy), 1.f, 500000.f);
}

void cb_key(GLFWwindow* w, int key, int, int act, int) {
    if (act != GLFW_PRESS && act != GLFW_REPEAT) return;
    if (key == GLFW_KEY_Q || key == GLFW_KEY_ESCAPE) glfwSetWindowShouldClose(w, 1);
    if (key == GLFW_KEY_R) g.cam = g.cam0;
}

void cb_resize(GLFWwindow*, int w, int h) {
    g.W = w; g.H = h; glViewport(0, 0, w, h);
}

// ── Main ──────────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {
    const char* file = argc > 1 ? argv[1] : "output.gcode";
    std::cout << "Loading: " << file << "\n";
    GcodeData gcode = load_gcode(file);
    std::cout << "Layers: " << gcode.layer_min << " – " << gcode.layer_max << "\n";
    std::cout << "Layer height: " << gcode.groups[0].tube_lh << " mm"
              << "  Line width: "  << gcode.groups[0].tube_lw << " mm\n";
    std::cout << "Retractions: " << gcode.retract_pts.size() << "\n";

    // Build retraction marker geometry (octahedra, 1.5× line width), grouped by (layer, seq)
    {
        float marker_r = gcode.groups[0].tube_lw * 1.5f;
        auto& rg = gcode.groups[ST_RETRACT];
        int   cur_layer = INT_MIN;
        int   cur_seq = -1;
        uint32_t lr_start = 0;
        for (auto& rpt : gcode.retract_pts) {
            if (rpt.layer != cur_layer || rpt.seq != cur_seq) {
                if (cur_layer != INT_MIN)
                    rg.layer_ranges.push_back({cur_layer, cur_seq, lr_start,
                                               (uint32_t)rg.idx.size() - lr_start});
                cur_layer = rpt.layer;
                cur_seq = rpt.seq;
                lr_start  = (uint32_t)rg.idx.size();
            }
            build_sphere_marker(rg.verts, rg.idx, rpt.pos, marker_r);
        }
        if (cur_layer != INT_MIN)
            rg.layer_ranges.push_back({cur_layer, cur_seq, lr_start,
                                       (uint32_t)rg.idx.size() - lr_start});
    }

    // Build per-layer sequence counts for the playback slider
    for (auto& cmd : gcode.draw_cmds)
        gcode.layer_seq_count[cmd.layer] = std::max(gcode.layer_seq_count[cmd.layer], cmd.seq + 1);

    if (!glfwInit()) { std::cerr << "glfwInit failed\n"; return 1; }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    GLFWwindow* win = glfwCreateWindow(g.W, g.H, "G-code Viewer", nullptr, nullptr);
    if (!win) { std::cerr << "Window creation failed\n"; glfwTerminate(); return 1; }
    glfwMakeContextCurrent(win);
    glfwSwapInterval(1);

    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) { std::cerr << "glewInit failed\n"; return 1; }

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // No primitive restart needed — tube meshes use plain GL_TRIANGLES

    GLuint prog     = make_program();
    GLint  mvp_loc  = glGetUniformLocation(prog, "MVP");
    GLint  col_loc  = glGetUniformLocation(prog, "color");
    GLint  ldir_loc = glGetUniformLocation(prog, "lightDir");

    // Upload tube geometry
    // Vertex layout: px py pz nx ny nz  (6 floats = sizeof(Vert))
    for (auto& grp : gcode.groups) {
        if (grp.verts.empty()) continue;
        glGenVertexArrays(1, &grp.vao);
        glGenBuffers(1, &grp.vbo);
        glGenBuffers(1, &grp.ebo);
        glBindVertexArray(grp.vao);
        glBindBuffer(GL_ARRAY_BUFFER, grp.vbo);
        glBufferData(GL_ARRAY_BUFFER,
            (GLsizeiptr)(grp.verts.size() * sizeof(Vert)),
            grp.verts.data(), GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, grp.ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,
            (GLsizeiptr)(grp.idx.size() * sizeof(uint32_t)),
            grp.idx.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vert), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vert), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(1);
        glBindVertexArray(0);

        std::cout << kNames[grp.type] << ": "
                  << grp.verts.size() << " verts, "
                  << grp.idx.size() / 3 << " triangles\n";
    }

    // Occluder at min radius
    float occ_r = (gcode.r_min < 1e8f ? gcode.r_min : 10.f) * 0.92f;
    Mesh  occ   = build_cylinder(occ_r, gcode.x_min, gcode.x_max);

    // Camera
    float xmid = (gcode.x_min + gcode.x_max) * 0.5f;
    float span  = std::max(gcode.x_max - gcode.x_min, gcode.r_max * 4.f);
    g.cam.target = glm::vec3(xmid, 0, 0);
    g.cam.dist   = span * 0.7f;
    g.cam0 = g.cam;

    for (int i = 0; i < ST_COUNT; ++i) g.visible[i] = true;
    g.total_layer_min = gcode.layer_min;
    g.total_layer_max = gcode.layer_max;
    g.show_layer_min  = gcode.layer_min;
    g.show_layer_max  = gcode.layer_max;

    glfwSetMouseButtonCallback(win, cb_mouse);
    glfwSetCursorPosCallback(win, cb_pos);
    glfwSetScrollCallback(win, cb_scroll);
    glfwSetKeyCallback(win, cb_key);
    glfwSetFramebufferSizeCallback(win, cb_resize);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(win, true);
    ImGui_ImplOpenGL3_Init("#version 330");
    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 6.f;
    style.FrameRounding  = 4.f;
    style.WindowPadding  = {12, 12};
    style.ItemSpacing    = {8, 6};

    // Fixed world-space light direction (upper-right-front)
    glm::vec3 ldir = glm::normalize(glm::vec3(0.4f, 1.0f, 0.6f));

    while (!glfwWindowShouldClose(win)) {
        glfwPollEvents();

        // ── ImGui frame ──
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::SetNextWindowPos({(float)g.W - 240.f, 8.f}, ImGuiCond_Always);
        ImGui::SetNextWindowSize({232.f, 0.f}, ImGuiCond_Always);
        ImGui::Begin("Controls", nullptr,
            ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_AlwaysAutoResize);

        ImGui::TextUnformatted("Layers");
        ImGui::SetNextItemWidth(-1);
        ImGui::DragIntRange2("##layers",
            &g.show_layer_min, &g.show_layer_max,
            1.f,
            g.total_layer_min, g.total_layer_max,
            "Min: %d", "Max: %d");

        // Single-layer row: [◄] [typed input] [►]  [All]
        float fh = ImGui::GetFrameHeight();
        float sp = ImGui::GetStyle().ItemSpacing.x;
        float all_w = ImGui::CalcTextSize("All").x + ImGui::GetStyle().FramePadding.x * 2;
        if (ImGui::ArrowButton("##lprev", ImGuiDir_Left)) {
            int v = std::max(g.total_layer_min, g.show_layer_min - 1);
            g.show_layer_min = g.show_layer_max = v;
        }
        ImGui::SameLine(0, 2);
        ImGui::SetNextItemWidth(-(fh + 2 + sp + all_w + sp));
        int snap = g.show_layer_min;
        if (ImGui::InputInt("##lsnap", &snap, 0, 0)) {
            snap = std::max(g.total_layer_min, std::min(g.total_layer_max, snap));
            g.show_layer_min = g.show_layer_max = snap;
        }
        ImGui::SameLine(0, 2);
        if (ImGui::ArrowButton("##lnext", ImGuiDir_Right)) {
            int v = std::min(g.total_layer_max, g.show_layer_max + 1);
            g.show_layer_min = g.show_layer_max = v;
        }
        ImGui::SameLine();
        if (ImGui::Button("All")) {
            g.show_layer_min = g.total_layer_min;
            g.show_layer_max = g.total_layer_max;
        }
        ImGui::TextDisabled("range %d – %d", g.total_layer_min, g.total_layer_max);

        // Playback slider — always shown; active only on a single layer
        ImGui::Separator();
        ImGui::TextUnformatted("Playback");
        bool single_layer = (g.show_layer_min == g.show_layer_max);
        if (single_layer) {
            int lyr = g.show_layer_min;
            int total_seq = gcode.layer_seq_count.count(lyr) ? gcode.layer_seq_count.at(lyr) : 0;
            if (lyr != g.seq_layer) {   // reset when layer changes
                g.seq_layer = lyr;
                g.show_seq  = total_seq - 1;
            }
            g.show_seq = std::max(0, std::min(g.show_seq, total_seq - 1));

            // Display path counter
            ImGui::Text("Path %d / %d", g.show_seq, total_seq);
            
            // Use float slider with power-law scaling for fine control on layers with many paths
            float slider_max = 1000.f;
            float power = 0.3f;  // < 1.0 gives more resolution at start

            // Convert show_seq to slider value
            float slider_val = total_seq > 1
                ? powf((float)g.show_seq / (float)(total_seq - 1), power) * slider_max
                : 0.f;
            
            float btn_w = ImGui::CalcTextSize("End").x + ImGui::GetStyle().FramePadding.x * 2;
            ImGui::SetNextItemWidth(-(btn_w + ImGui::GetStyle().ItemSpacing.x));
            
            if (ImGui::SliderFloat("##seq", &slider_val, 0.f, slider_max, "", 
                                   ImGuiSliderFlags_NoInput)) {
                // Convert slider value back to show_seq
                float t = slider_val / slider_max;
                g.show_seq = (int)(powf(t, 1.f / power) * (float)(total_seq - 1) + 0.5f);
                g.show_seq = std::max(0, std::min(g.show_seq, total_seq - 1));
            }
            ImGui::SameLine();
            if (ImGui::Button("End")) {
                g.show_seq = total_seq - 1;
            }
        } else {
            g.seq_layer = INT_MIN;
            g.show_seq  = -1;
            ImGui::TextDisabled("(select a single layer)");
        }

        ImGui::Separator();
        ImGui::TextUnformatted("Line types");

        bool all_on = true;
        for (int i = 0; i < ST_COUNT; ++i) if (!g.visible[i]) { all_on = false; break; }
        if (ImGui::Checkbox("All##all", &all_on))
            for (int i = 0; i < ST_COUNT; ++i) g.visible[i] = all_on;

        for (int i = 0; i < ST_COUNT; ++i) {
            ImGui::ColorButton(("##col" + std::to_string(i)).c_str(),
                ImVec4(kColors[i][0], kColors[i][1], kColors[i][2], 1.f),
                ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoBorder,
                ImVec2(14, 14));
            ImGui::SameLine();
            ImGui::Checkbox(kNames[i], &g.visible[i]);
        }

        ImGui::Separator();
        ImGui::TextDisabled("L-drag: orbit");
        ImGui::TextDisabled("R-drag: pan");
        ImGui::TextDisabled("Scroll: zoom");
        ImGui::TextDisabled("R: reset camera");

        ImGui::End();

        // ── G-code feed window (single-layer playback only) ──
        bool seq_feed_mode = (g.show_layer_min == g.show_layer_max && g.show_seq >= 0);
        if (seq_feed_mode) {
            float feed_w = 420.f;
            float feed_h = (float)g.H * 0.30f;
            ImGui::SetNextWindowPos({8.f, (float)g.H - feed_h - 8.f}, ImGuiCond_Always);
            ImGui::SetNextWindowSize({feed_w, feed_h}, ImGuiCond_Always);
            ImGui::Begin("G-code Feed", nullptr,
                ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize |
                ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_HorizontalScrollbar);

            bool seq_changed = (g.show_seq != g.feed_prev_seq);
            g.feed_prev_seq = g.show_seq;

            int lyr = g.show_layer_min;
            for (auto& cmd : gcode.draw_cmds) {
                if (cmd.layer != lyr) continue;
                if (cmd.seq > g.show_seq) break;
                bool is_current = (cmd.seq == g.show_seq);
                const float* col = kColors[cmd.type];
                for (auto& ln : cmd.lines) {
                    if (is_current)
                        ImGui::TextColored(ImVec4(col[0], col[1], col[2], 1.f), "%s", ln.c_str());
                    else
                        ImGui::TextDisabled("%s", ln.c_str());
                }
            }

            if (seq_changed)
                ImGui::SetScrollHereY(1.0f);

            ImGui::End();
        } else {
            g.feed_prev_seq = -2;
        }

        // ── 3D render ──
        glClearColor(0.12f, 0.12f, 0.12f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glm::mat4 proj = glm::perspective(glm::radians(45.f), (float)g.W / (float)g.H, 0.1f, 100000.f);
        glm::mat4 mvp  = proj * g.cam.view();

        glUseProgram(prog);
        glUniformMatrix4fv(mvp_loc, 1, GL_FALSE, glm::value_ptr(mvp));
        glUniform3fv(ldir_loc, 1, glm::value_ptr(ldir));

        // Occluder
        glDepthMask(GL_TRUE);
        glUniform4f(col_loc, 0.f, 0.f, 0.f, 1.f);
        glBindVertexArray(occ.vao);
        glDrawElements(GL_TRIANGLE_STRIP, occ.count, GL_UNSIGNED_INT, nullptr);
        glBindVertexArray(0);

        // Opaque passes first, then transparent (travel etc.)
        bool seq_mode = (g.show_layer_min == g.show_layer_max && g.show_seq >= 0);
        for (int pass = 0; pass < 2; ++pass) {
            if (seq_mode) {
                // Sequential playback: iterate draw_cmds in order, stop at show_seq
                GLuint bound_vao = 0;
                int    bound_i   = -1;
                for (auto& cmd : gcode.draw_cmds) {
                    if (cmd.layer != g.show_layer_min) continue;
                    if (cmd.seq > g.show_seq) continue;
                    if (!g.visible[cmd.type]) continue;
                    bool transparent = (kAlpha[cmd.type] < 1.0f);
                    if (pass == 0 &&  transparent) continue;
                    if (pass == 1 && !transparent) continue;
                    auto& grp = gcode.groups[cmd.type];
                    if (!grp.vao) continue;
                    if (pass == 1) glDepthMask(GL_FALSE);
                    if ((int)cmd.type != bound_i) {
                        glUniform4f(col_loc, kColors[cmd.type][0], kColors[cmd.type][1],
                                    kColors[cmd.type][2], kAlpha[cmd.type]);
                        bound_i = (int)cmd.type;
                    }
                    if (grp.vao != bound_vao) { glBindVertexArray(grp.vao); bound_vao = grp.vao; }
                    glDrawElements(GL_TRIANGLES, (GLsizei)cmd.idx_count, GL_UNSIGNED_INT,
                                   (void*)(uintptr_t)(cmd.idx_offset * sizeof(uint32_t)));
                    if (pass == 1) glDepthMask(GL_TRUE);
                }
                if (bound_vao) glBindVertexArray(0);

                // Retraction markers for this layer (not part of draw_cmds)
                if (pass == 0 && g.visible[ST_RETRACT]) {
                    auto& rg = gcode.groups[ST_RETRACT];
                    if (rg.vao) {
                        glUniform4f(col_loc, kColors[ST_RETRACT][0], kColors[ST_RETRACT][1],
                                    kColors[ST_RETRACT][2], kAlpha[ST_RETRACT]);
                        glBindVertexArray(rg.vao);
                        for (auto& lr : rg.layer_ranges) {
                            if (lr.layer_num != g.show_layer_min) continue;
                            // In sequential mode, only show retractions up to current sequence
                            if (seq_mode && lr.seq > g.show_seq) continue;
                            glDrawElements(GL_TRIANGLES, (GLsizei)lr.idx_count, GL_UNSIGNED_INT,
                                           (void*)(uintptr_t)(lr.idx_offset * sizeof(uint32_t)));
                        }
                        glBindVertexArray(0);
                    }
                }
            } else {
                // Normal mode: draw by layer range
                for (int i = 0; i < ST_COUNT; ++i) {
                    bool transparent = (kAlpha[i] < 1.0f);
                    if (pass == 0 &&  transparent) continue;
                    if (pass == 1 && !transparent) continue;
                    if (!g.visible[i]) continue;
                    auto& grp = gcode.groups[i];
                    if (!grp.vao) continue;
                    if (transparent) glDepthMask(GL_FALSE);
                    glUniform4f(col_loc, kColors[i][0], kColors[i][1], kColors[i][2], kAlpha[i]);
                    glBindVertexArray(grp.vao);
                    for (auto& lr : grp.layer_ranges) {
                        if (lr.layer_num < g.show_layer_min || lr.layer_num > g.show_layer_max) continue;
                        glDrawElements(GL_TRIANGLES, (GLsizei)lr.idx_count, GL_UNSIGNED_INT,
                                       (void*)(uintptr_t)(lr.idx_offset * sizeof(uint32_t)));
                    }
                    glBindVertexArray(0);
                    if (transparent) glDepthMask(GL_TRUE);
                }
            }
        }

        // Draw playback position pointer in sequential mode
        if (seq_mode && g.visible[ST_RETRACT]) {
            // Find the last drawn command to get the current position
            Vec3 current_pos = {0, 0, 0};
            bool found = false;
            for (auto it = gcode.draw_cmds.rbegin(); it != gcode.draw_cmds.rend(); ++it) {
                if (it->layer == g.show_layer_min && it->seq <= g.show_seq && g.visible[it->type]) {
                    current_pos = it->end_pos;
                    found = true;
                    break;
                }
            }
            if (found) {
                // Build and draw cone marker at current position
                std::vector<Vert> cone_verts;
                std::vector<uint32_t> cone_idx;
                float cone_len = 3.f;  // mm
                float cone_rad = 1.5f; // mm
                // Point in +X direction (axial)
                build_cone_marker(cone_verts, cone_idx, current_pos, Vec3{1, 0, 0}, cone_len, cone_rad);
                
                // Upload and draw
                GLuint cone_vao, cone_vbo, cone_ebo;
                glGenVertexArrays(1, &cone_vao);
                glGenBuffers(1, &cone_vbo);
                glGenBuffers(1, &cone_ebo);
                glBindVertexArray(cone_vao);
                glBindBuffer(GL_ARRAY_BUFFER, cone_vbo);
                glBufferData(GL_ARRAY_BUFFER,
                    (GLsizeiptr)(cone_verts.size() * sizeof(Vert)),
                    cone_verts.data(), GL_DYNAMIC_DRAW);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cone_ebo);
                glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                    (GLsizeiptr)(cone_idx.size() * sizeof(uint32_t)),
                    cone_idx.data(), GL_DYNAMIC_DRAW);
                glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vert), (void*)0);
                glEnableVertexAttribArray(0);
                glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vert), (void*)(3 * sizeof(float)));
                glEnableVertexAttribArray(1);
                
                // Draw with bright white color, fully opaque
                glDepthMask(GL_TRUE);
                glUniform4f(col_loc, 1.f, 1.f, 0.f, 1.f);  // Yellow
                glBindVertexArray(cone_vao);
                glDrawElements(GL_TRIANGLES, (GLsizei)cone_idx.size(), GL_UNSIGNED_INT, nullptr);
                glBindVertexArray(0);
                
                glDeleteVertexArrays(1, &cone_vao);
                glDeleteBuffers(1, &cone_vbo);
                glDeleteBuffers(1, &cone_ebo);
            }
        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(win);
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(win);
    glfwTerminate();
    return 0;
}
