function [ fun ] = geom_cyl_xy( data, inside )
%GEOM_CYL Summary of this function goes here
%   Detailed explanation goes here
% data(1) is x center
%data(2) is y center
%data(3) is R
%hole_flag==1 then the outer of cylinder is taken

fun=@(xyz) sum( ~xor( (xyz(1)-data(1)).^2 + (xyz(2)-data(2)).^2 <= data(3).^2   ,inside));
end

