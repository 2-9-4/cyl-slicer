module subdivided_cube(size=1, center=false, dim=5) {
    // Standardize size to a 3-element vector [x, y, z]
    s = is_list(size) ? size : [size, size, size];
    
    // Calculate the translation offset if center is true
    offset = center ? [-s[0]/2, -s[1]/2, -s[2]/2] : [0, 0, 0];
    
    step_x = s[0] / dim;
    step_y = s[1] / dim;
    step_z = s[2] / dim;

    translate(offset)
    for (x = [0 : step_x : s[0] - 0.001]) {
        for (y = [0 : step_y : s[1] - 0.001]) {
            for (z = [0 : step_z : s[2] - 0.001]) {
                // 0.01 overlap ensures a single manifold object
                translate([x, y, z]) 
                    cube([step_x + 0.01, step_y + 0.01, step_z + 0.01]);
            }
        }
    }
}

translate([60+35,0,0])rotate([35,0,45])subdivided_cube(20,center=true);