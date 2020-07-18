function [ out ] = bin2xy( b,params)
%INDICES_TRANSLATION INTO XY
%decypher 1-by-(dz+dy) binary line 
%into 1_by_2 decimal array of 2D spatial coordinates
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
hx=params.hx;
hy=params.hy;
  
y=bi2de(b(:,1:dy));
x=bi2de(b(:,dy+1:dy+dx));
out=[(x)*hx (y)*hy];

end

