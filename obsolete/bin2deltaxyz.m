function [ out ] = bin2deltaxyz( b,params)
%INDICES_TRANSLATION INTO DELTAS
%decypher 1-by-(dz+dy+2) binary line 
%into 1_by_2 decimal array of 2D spatial coordinates in Deltas space
% equispaced square structured lattice implied

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
  
z=bi2de(b(:,1:dz+1));
y=bi2de(b(:,dz+2:dz+dy+2));
x=bi2de(b(:,dz+dy+3:dz+dy+dx));

out=[(x-2^dx)*hx (y-2^dy)*hy (z-2^dz)*hz];

end

