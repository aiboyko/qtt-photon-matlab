function mask_of_coords = geom_box(params,ix1,ix2,iy1,iy2,iz1,iz2)
%Lossless version of box-generation
% works in index-numbers space

%not sure if indexing is entirely correct
 dx=params.grid.dx;
 dy=params.grid.dy;
 dz=params.grid.dz;
 hx=params.grid.hx;
 hy=params.grid.hy;
 hz=params.grid.hz;
 tol=params.tol;

if ix1<0
   ix1=1;
end
tth1x=tt_heaviside(2,dx,ix1);
if ix2>=(2^dx)
    tth2x=tt_ones(2,dx);
else
    tth2x=1-tt_heaviside(2,dx,ix2+1);
end
tth3x=round(tth1x.*tth2x,tol);

if iy1<0
   iy1=1;
end
tth1y=tt_heaviside(2,dy,iy1);
if iy2>=(2^dy)
    tth2y=tt_ones(2,dy);
else
    tth2y=1-tt_heaviside(2,dy,iy2+1);
end
tth3y=round(tth1y.*tth2y,tol);

if iz1<0
    iz1=1;
end

if dz>0
    tth1z=tt_heaviside(2,dz,iz1);
    if iz2>=(2^dz)
        tth2z=tt_ones(2,dz);
    else
        tth2z=1-tt_heaviside(2,dz,iz2+1); %then iz2 is still 1
    end
    tth3z=round(tth1z.*tth2z,tol);
else
    tth3z=tt_tensor(1);   
end
mask_of_coords=round(purify(mtkron(tth3z,tth3y,tth3x)),tol);
end