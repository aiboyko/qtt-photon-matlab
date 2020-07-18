function [ v_reshaped,X,Y,Z ] = F1_LtoF1_3D( v,params )
%QTT2FIELD Summary of this function goes here
dx=params.dx;
dy=params.dy;
dz=params.dz;
hx=params.hx;
hy=params.hy;
hz=params.hz;
tol=params.tol;
 
v_reshaped=reshape(v,[2^dz 2^dy 2^dx]);
v_reshaped=permute(v_reshaped,[3 2 1]);

[Y,X,Z]=meshgrid(hy*(0:(2^dy-1)),hx*(0:(2^dx-1)),hz*(0:(2^dz-1)));

end

