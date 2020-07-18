function [ out ] = bin2xyz( b,params,field )
%INDICES_TRANSLATION 
%decypher M_by_(dz+dy+dz) binary matrix
%into M_by_3 decimal array of 3D spatial coordinates
% equispaced structured lattice implied

%it looks like it works

%#0 [ 0 0 0 .. 0 0 0 .. 1 0 1]
%       ^   ..   ^   ..   ^
%       dz       dy       dx
%#1 [ 1 0 .. 0]
%..
%#M [ 1 1 1 .. 1 1 1 .. 1 1 1]

dx=params.dx;
dy=params.dy;
dz=params.dz;
hx=params.hx;
hy=params.hy;
hz=params.hz;

   if nargin>7
       delta=deltaXYZ(field);
   else
       delta=[ 0 0 0];
   end
   
z=bi2de(b(:,1:dz))*hz;

y=bi2de(b(:,dz+1:dz+dy))*hy;

x=bi2de(b(:,dz+dy+1:dz+dy+dx))*hx;

out=[x y z]+delta.*[hx,hy,hz];

end

