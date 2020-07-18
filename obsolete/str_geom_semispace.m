function [ strfun ] = str_geom_semispace( xc ,n )

%GEOM_CYL Summary of this function goes here
%   Detailed explanation goes here
%n is a column vector
%hole_flag==1 then the outer of cylinder is taken

strfun=strcat(num2str(n(1)),'*(xyz(1)-',num2str(xc(1)),')+',num2str(n(2)),...
    '*(xyz(2)-',num2str(xc(2)),')+',num2str(n(3)),'*(xyz(3)-',num2str(xc(3)),')>=0');

end



