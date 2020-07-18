function [ fun ] = geom_semispace( xc ,n )

%GEOM_CYL Summary of this function goes here
%   Detailed explanation goes here
%n is a column vector
%hole_flag==1 then the outer of cylinder is taken

fun=@(xyz)(xyz-xc)*(n.') >=0;

end



