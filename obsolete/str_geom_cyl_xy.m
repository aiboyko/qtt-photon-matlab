function [ strfun ] = str_geom_cyl_xy( data, inside )
%GEOM_CYL Summary of this function goes here
%   Detailed explanation goes here
% data(1) is x center
%data(2) is y center
%data(3) is R
%hole_flag==1 then the outer of cylinder is taken

strfun=strcat('~xor((xyz(1)-',num2str(data(1)),').^2 + (xyz(2)-',...
    num2str(data(2)),').^2 <=(',num2str(data(3)),').^2,',num2str(inside),')');
end

