function [ bawks ] = geom_box( x1, x2 ,y1 ,y2, z1, z2 )

%GEOM_CYL Summary of this function goes here
%   Detailed explanation goes here
%data(1) is x center
%data(2) is y center
%data(3) is R
%hole_flag==1 then the outer of cylinder is taken

p1=geom_semispace([x1 0 0],[1 0 0]);
p2=geom_semispace([x2 0 0],[-1 0 0]);
p3=geom_semispace([0 y1  0],[0 1 0]);
p4=geom_semispace([0 y2 0],[0 -1 0]);
p5=geom_semispace([0 0 z1],[0 0 1]);
p6=geom_semispace([0 0 z2],[0 0 -1]);


bawks=@(xyz)p1(xyz)&&p2(xyz)&&p3(xyz)&&p4(xyz)&&p5(xyz)&&p6(xyz);
end





