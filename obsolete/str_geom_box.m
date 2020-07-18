function [ strfun ] = str_geom_box( x1, x2 ,y1 ,y2, z1, z2 )

%GEOM_CYL Summary of this function goes here
%   Detailed explanation goes here
%data(1) is x center
%data(2) is y center
%data(3) is R
%hole_flag==1 then the outer of cylinder is taken

C={...
str_geom_semispace([x1 0 0],[1  0  0]),...
str_geom_semispace([x2 0 0],[-1 0  0]),...
str_geom_semispace([0 y1 0],[0  1  0]),...
str_geom_semispace([0 y2 0],[0 -1  0]),...
str_geom_semispace([0 0 z1],[0  0  1]),...
str_geom_semispace([0 0 z2],[0  0 -1])};



strfun=strjoin(C,'&&');
end





