// Cylindrical retraction coupon. Mesh Z is axial; XY is the radial plane.
// Slicer output Z is radial height above bed_radius, NOT this mesh's Z.
// openscad -o /tmp/retraction_posts.stl calibration/retraction_posts.scad
bed_radius = 30;
base_height = 2;
band_height = 3;
bands = 5;
post_width = 6;              // axial and arc width at the build surface (mm)
axial_spacing = 22;          // centre to centre (mm)
arc_spacing = 20;            // centre to centre at bed_radius (mm)
margin = 3;
layout = "both";            // "axial", "rotary", or "both" (four posts)
// Optional witness grooves. Default off to avoid changing geometry between bands.
// If enabled, do not score the transition shells around grooves.
groove_depth = 0;
groove_height = 0.4;
$fn = 720;

assert(bed_radius > 0 && base_height > 0 && band_height > groove_height);
assert(bands >= 1 && bands == floor(bands));
assert(post_width > 2*groove_depth && axial_spacing > post_width && arc_spacing > post_width);
assert(layout == "axial" || layout == "rotary" || layout == "both");

// rotate_extrude revolves an (r, axial) rectangle around mesh Z.
module sector(r0, r1, theta, arc_width, axial, length) {
    angle = arc_width / bed_radius * 180 / PI;
    assert(angle > 0 && angle < 360);
    rotate([0, 0, theta - angle/2])
        rotate_extrude(angle=angle)
            translate([r0, axial-length/2]) square([r1-r0, length]);
}

total_height = base_height + bands * band_height;
axials = layout == "rotary" ? [0] : [-axial_spacing/2, axial_spacing/2];
arcs = layout == "axial" ? [0] : [-arc_spacing/2, arc_spacing/2];
base_arc = (layout == "axial" ? 0 : arc_spacing) + post_width + 2*margin;
base_length = (layout == "rotary" ? 0 : axial_spacing) + post_width + 2*margin;

union() {
    sector(bed_radius, bed_radius + base_height, 90, base_arc, 0, base_length);
    for (axial = axials, arc = arcs) {
        theta = 90 + arc / bed_radius * 180 / PI;
        // Core passes through all grooves, and overlaps the curved base.
        sector(bed_radius + base_height - 0.1, bed_radius + total_height,
               theta, post_width-2*groove_depth, axial, post_width-2*groove_depth);
        for (band = [0:bands-1]) {
            r0 = bed_radius + base_height + band*band_height;
            r1 = r0 + band_height - (band < bands-1 ? groove_height : 0);
            sector(r0, r1, theta, post_width, axial, post_width);
        }
    }
}
