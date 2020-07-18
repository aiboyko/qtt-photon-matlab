function [ idxx ] = xyz2bin( xyz,params,field )
%INDICES_TRANSLATION 
%encode M_by_3 decimal array of 3D spatial coordinates
%into M_by_(dx+dy+dz) binary matrix
%into M_by_3 decimal array of 3D spatial coordinates
% structured lattice implied

%#0 [ 0 0 0 .. 0 0 0 .. 1 0 1]
%       ^   ..   ^   ..   ^
%       dz       dy       dx
%#1 [ 1 0 .. 0]
%..
%#M [ 1 1 1 .. 1 1 1 .. 1 1 1]

%it looks like it works


dx=params.dx;
dy=params.dy;
dz=params.dz;
hx=params.hx;
hy=params.hy;
hz=params.hz;

   if nargin>2
       delta=deltaXYZ(field);
   else
       delta=[ 0 0 0];
   end
   center_xyz=floor((xyz-delta)./[ hx hy hz]);

if dz>0
    idxx=[de2bi(center_xyz(3),dz) ,de2bi(center_xyz(2),dy) ,...
        de2bi(center_xyz(1),dx)];
else
    idxx=[de2bi(center_xyz(2),dy) ,...
    de2bi(center_xyz(1),dx)];
end  
%this is a purely binary(0 or 1) encoder. To use in ttv one need to add 1.

end

