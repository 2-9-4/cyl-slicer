// gcode_viewer.cpp — High-performance OpenGL G-code viewer for cylindrical printer
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
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// ── Section types ─────────────────────────────────────────────────────────────

enum SType {
    ST_EXTERNAL = 0, ST_INTERNAL, ST_SOLID_INFILL,
    ST_RADIAL_FULL,  ST_RADIAL,   ST_AXIAL,
    ST_SUPPORT,      ST_TRAVEL,   ST_UNKNOWN, ST_COUNT
};

static const char* kNames[ST_COUNT] = {
    "External Perimeter", "Internal Perimeter", "Solid Infill",
    "Infill Radial Full", "Infill Radial",      "Infill Axial",
    "Support",            "Travel",             "Unknown",
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
    {0.533f, 0.533f, 0.533f},  // UNKNOWN       #888888
};

static const float kAlpha[ST_COUNT] = {
    1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.10f, 0.6f
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

// ── Geometry ──────────────────────────────────────────────────────────────────

struct Vec3 { float x, y, z; };

Vec3 cyl_to_cart(float theta_deg, float ax, float r) {
    float t = std::fmod(theta_deg, 360.f) * (float)M_PI / 180.f;
    return { ax, r * std::cos(t), r * std::sin(t) };
}

struct LayerRange {
    int      layer_num;
    uint32_t idx_offset;  // element offset into index buffer
    uint32_t idx_count;
};

struct PathGroup {
    SType                 type;
    std::vector<Vec3>     verts;
    std::vector<uint32_t> idx;          // 0xFFFFFFFF = primitive restart
    std::vector<LayerRange> layer_ranges;
    GLuint vao = 0, vbo = 0, ebo = 0;
};

struct GcodeData {
    PathGroup groups[ST_COUNT];
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

    float cur_x = 0.f, cur_t = 0.f;
    float radius = 10.f;
    bool  has_pos = false;
    bool  in_ext  = false;
    SType ext_type = ST_UNKNOWN;
    SType sec_type = ST_UNKNOWN;
    int   cur_layer = 0;

    uint32_t layer_idx_start[ST_COUNT] = {};

    auto push_vert = [](PathGroup& g, Vec3 v) {
        uint32_t i = (uint32_t)g.verts.size();
        g.verts.push_back(v);
        g.idx.push_back(i);
    };
    auto start_path = [](PathGroup& g) {
        if (!g.idx.empty()) g.idx.push_back(0xFFFFFFFF);
    };
    auto close_layer = [&]() {
        for (int i = 0; i < ST_COUNT; ++i) {
            uint32_t start = layer_idx_start[i];
            uint32_t end   = (uint32_t)d.groups[i].idx.size();
            if (end > start)
                d.groups[i].layer_ranges.push_back({cur_layer, start, end - start});
            layer_idx_start[i] = end;
        }
    };

    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        if (line[0] == ';') {
            if (line.rfind(";LAYER:", 0) == 0) {
                close_layer();
                std::string layer_str = line.substr(7);
                // Accept integers (positive or negative); skip non-numeric identifiers
                if (!layer_str.empty() && (std::isdigit(layer_str[0]) || layer_str[0] == '-')) {
                    cur_layer = std::stoi(layer_str);
                    d.layer_min = std::min(d.layer_min, cur_layer);
                    d.layer_max = std::max(d.layer_max, cur_layer);
                }
                in_ext = false;
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
        float nx = cur_x, nt = cur_t;

        std::string tok;
        while (ss >> tok) {
            if (tok.size() < 2) continue;
            try {
                float v = std::stof(tok.substr(1));
                if      (tok[0] == 'X') nx = v;
                else if (tok[0] == 'Y') nt = v;
            } catch (...) {}
        }

        d.x_min = std::min(d.x_min, nx);
        d.x_max = std::max(d.x_max, nx);
        Vec3 nv = cyl_to_cart(nt, nx, radius);

        if (rapid) {
            in_ext = false;
            if (has_pos) {
                auto& g = d.groups[ST_TRAVEL];
                start_path(g);
                push_vert(g, cyl_to_cart(cur_t, cur_x, radius));
                push_vert(g, nv);
            }
        } else {
            if (!in_ext || sec_type != ext_type) {
                in_ext   = true;
                ext_type = sec_type;
                auto& g = d.groups[ext_type];
                start_path(g);
                if (has_pos) push_vert(g, cyl_to_cart(cur_t, cur_x, radius));
            }
            push_vert(d.groups[ext_type], nv);
        }

        cur_x = nx; cur_t = nt; has_pos = true;
    }

    close_layer();  // flush last layer
    if (d.layer_min == INT_MAX) { d.layer_min = 0; d.layer_max = 0; }
    return d;
}

// ── Cylinder mesh (occluder) ──────────────────────────────────────────────────

struct Mesh { GLuint vao = 0, vbo = 0, ebo = 0; int count = 0; };

Mesh build_cylinder(float r, float x0, float x1, int slices = 64) {
    std::vector<float>    v;
    std::vector<uint32_t> idx;

    for (int xi = 0; xi < 2; ++xi) {
        float ax = xi == 0 ? x0 : x1;
        for (int i = 0; i <= slices; ++i) {
            float t = 2.f * (float)M_PI * i / slices;
            v.push_back(ax);
            v.push_back(r * std::cos(t));
            v.push_back(r * std::sin(t));
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
    glBufferData(GL_ARRAY_BUFFER, v.size() * sizeof(float), v.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, idx.size() * sizeof(uint32_t), idx.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 12, nullptr);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
    m.count = (int)idx.size();
    return m;
}

// ── Shaders ───────────────────────────────────────────────────────────────────

static const char* kVert = R"glsl(
#version 330 core
layout(location=0) in vec3 pos;
uniform mat4 MVP;
void main() { gl_Position = MVP * vec4(pos, 1.0); }
)glsl";

static const char* kFrag = R"glsl(
#version 330 core
uniform vec4 color;
out vec4 fragColor;
void main() { fragColor = color; }
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
    glEnable(GL_PRIMITIVE_RESTART);
    glPrimitiveRestartIndex(0xFFFFFFFF);
    glLineWidth(2.0f);

    GLuint prog    = make_program();
    GLint  mvp_loc = glGetUniformLocation(prog, "MVP");
    GLint  col_loc = glGetUniformLocation(prog, "color");

    // Upload path geometry
    for (auto& grp : gcode.groups) {
        if (grp.verts.empty()) continue;
        glGenVertexArrays(1, &grp.vao);
        glGenBuffers(1, &grp.vbo);
        glGenBuffers(1, &grp.ebo);
        glBindVertexArray(grp.vao);
        glBindBuffer(GL_ARRAY_BUFFER, grp.vbo);
        glBufferData(GL_ARRAY_BUFFER, grp.verts.size() * sizeof(Vec3), grp.verts.data(), GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, grp.ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, grp.idx.size() * sizeof(uint32_t), grp.idx.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vec3), nullptr);
        glEnableVertexAttribArray(0);
        glBindVertexArray(0);
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

    // Register callbacks BEFORE ImGui init so ImGui chains onto them
    glfwSetMouseButtonCallback(win, cb_mouse);
    glfwSetCursorPosCallback(win, cb_pos);
    glfwSetScrollCallback(win, cb_scroll);
    glfwSetKeyCallback(win, cb_key);
    glfwSetFramebufferSizeCallback(win, cb_resize);

    // ImGui init — installs its own callbacks and chains to the ones above
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

        // Layer range slider
        ImGui::TextUnformatted("Layers");
        ImGui::SetNextItemWidth(-1);
        ImGui::DragIntRange2("##layers",
            &g.show_layer_min, &g.show_layer_max,
            1.f,
            g.total_layer_min, g.total_layer_max,
            "Min: %d", "Max: %d");
        // Quick "all layers" button
        if (ImGui::SmallButton("All")) {
            g.show_layer_min = g.total_layer_min;
            g.show_layer_max = g.total_layer_max;
        }
        ImGui::SameLine();
        ImGui::TextDisabled("(%d – %d)", g.total_layer_min, g.total_layer_max);

        ImGui::Separator();
        ImGui::TextUnformatted("Line types");

        // Checkboxes with color swatches
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

        // ── 3D render ──
        glClearColor(0.12f, 0.12f, 0.12f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glm::mat4 proj = glm::perspective(glm::radians(45.f), (float)g.W / (float)g.H, 0.1f, 100000.f);
        glm::mat4 mvp  = proj * g.cam.view();

        glUseProgram(prog);
        glUniformMatrix4fv(mvp_loc, 1, GL_FALSE, glm::value_ptr(mvp));

        // Occluder
        glDepthMask(GL_TRUE);
        glUniform4f(col_loc, 1.f, 1.f, 1.f, 1.f);
        glBindVertexArray(occ.vao);
        glDrawElements(GL_TRIANGLE_STRIP, occ.count, GL_UNSIGNED_INT, nullptr);
        glBindVertexArray(0);

        // Paths — draw only visible layers
        for (int i = 0; i < ST_COUNT; ++i) {
            if (!g.visible[i]) continue;
            auto& grp = gcode.groups[i];
            if (!grp.vao) continue;
            glUniform4f(col_loc, kColors[i][0], kColors[i][1], kColors[i][2], kAlpha[i]);
            glBindVertexArray(grp.vao);
            for (auto& lr : grp.layer_ranges) {
                if (lr.layer_num < g.show_layer_min || lr.layer_num > g.show_layer_max) continue;
                glDrawElements(GL_LINE_STRIP, (GLsizei)lr.idx_count, GL_UNSIGNED_INT,
                               (void*)(uintptr_t)(lr.idx_offset * sizeof(uint32_t)));
            }
            glBindVertexArray(0);
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
