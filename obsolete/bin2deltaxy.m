function [ out ] = bin2deltaxy( b,params)
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

%this function is created to explore G(k0|r1-r2|) over irreducible regular space of deltas, since it
%is exactly enough for qtt_bttb matrix constructor to create an reducible G matrix

dx=params.dx;
dy=params.dy;
hx=params.hx;
hy=params.hy;
  
y_idx=bi2de(b(:,1:dy+1));
x_idx=bi2de(b(:,dy+2:dy+dx+2));

out=[(x_idx-2^dx)*hx (y_idx-2^dy)*hy];
end

