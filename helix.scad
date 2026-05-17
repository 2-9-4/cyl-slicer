module tracePath(path, size){

for(i=[0:len(path)-1]){

    if(i+1 < len(path)){

        hull(){

            translate(path[i])

            rotate([0,0,atan2(path[i][1],path[i][0])])cube(size);

            translate(path[i+1])

            rotate([0,0,atan2(path[i+1][1],path[i+1][0])])cube(size);

        }

    }

}

}


function spiral(twist,height,spiral_radius) = [ for(i=[0:5:twist]) [spiral_radius*sin(i),spiral_radius*cos(i),i/twist*height]];

path = spiral(1200*2,50,35);

tracePath(path, 5);