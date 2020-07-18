function f=pointsource(x,y,z,params)
%this thing generates 2D hardcoded Bessel source
 dx=params.dx;
 dy=params.dy;
 dz=params.dz;
 hx=params.hx;
 hy=params.hy;
 hz=params.hz;
 tol=params.tol;
 k0=params.k0;
r=sqrt((x-0.5-.7*2^(dx-1)).^2+(y-0.5-.5*2^(dy-1)).^2);
f= 1/(r+1e-10).*exp(1j*k0*r);
end