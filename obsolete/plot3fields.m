function fig=plot3fields(F3_L, bf3_L, fullmask_3D,zcoord, params,nfig,name)
%this function plots 3 fields in xy-plane, that has a coordinate zcoord

 dx=params.dx;
 dy=params.dy;
 dz=params.dz;
%  hx=params.hx;
%  hy=params.hy;
%  hz=params.hz;
%  tol=params.tol;
%  k0=params.k0;
%  
N=2^(dx+dy+dz);
Fx_s_3d = F1_LtoF1_3D(F3_L(    1  :  N  ),params);
Fy_s_3d = F1_LtoF1_3D(F3_L( (N+1) :(2*N)),params);
Fz_s_3d = F1_LtoF1_3D(F3_L((2*N+1):(3*N)),params);
Jx_3d   = F1_LtoF1_3D(bf3_L(    (1):(N))  ,params);
Jy_3d   = F1_LtoF1_3D(bf3_L(  (N+1):(2*N)),params);
Jz_3d   = F1_LtoF1_3D(bf3_L((2*N+1):(3*N)),params);

if nargin>=11
    fig=figure(nfig);
else
    fig=figure();
end
if nargin==12
fig.Name=name;
end

subplot(2,3,1);
imagesc(real(Fx_s_3d(:,:,zcoord)));
title('Ex')
colorbar;
axis equal tight
hold on
contour(fullmask_3D(:,:,zcoord),1,'color','Black');
hold off

subplot(2,3,2);
imagesc(real(Fy_s_3d(:,:,zcoord)));
title('Ey')
colorbar;
axis equal tight
hold on
contour(fullmask_3D(:,:,zcoord),1,'color','Black');
hold off

subplot(2,3,3);
imagesc(real(Fz_s_3d(:,:,zcoord)));
title('Ez')
colorbar;
axis equal tight
hold on
contour(fullmask_3D(:,:,zcoord),1,'color','Black');
hold off

subplot(2,3,4);
imagesc(imag(Jx_3d(:,:,zcoord)));
colorbar;
axis equal tight
title('Jx')
hold on
contour(fullmask_3D(:,:,zcoord),1,'color','Black');
hold off

subplot(2,3,5);
imagesc(imag(Jy_3d(:,:,zcoord)));
colorbar;
title('Jy')
axis equal tight
hold on
contour(fullmask_3D(:,:,zcoord),1,'color','Black');
hold off

subplot(2,3,6);
imagesc(imag(Jz_3d(:,:,zcoord)));
colorbar;
title('Jz')
axis equal tight
hold on
contour(fullmask_3D(:,:,zcoord),1,'color','Black');
hold off
% % 
colormap(jet);