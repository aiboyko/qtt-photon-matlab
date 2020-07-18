%PEC
%DOESNT WORK!!
dx_mie=6;
dy_mie=dx_mie;
Lx=1;Ly=1;
hx=Lx/2^dx_mie;
hy=Ly/2^dy_mie;
EF=0.25;
alpha=0;%pi/2;

eps=3;
k0=2*pi/EF;
k1=k0;
k2=sqrt(eps)*k0;

[XX,YY]=meshgrid( hx/2+(0:2^dx_mie-1)*hx,hy/2+(0:2^dy_mie-1)*hy);
Xcenter=(Lx)/2; %[ . ][ . ][ . ][ . ] or  ._._._.
Ycenter=(Lx)/2;
XX=XX-Xcenter;
YY=YY-Ycenter;

[Th,R]=cart2pol(XX,YY);

field_0=exp(1j*k0*cos(alpha+Th+pi/2).*R);
% figure;imagesc(real(field_0));axis equal tight;  colormap(jet(2^13));
% 
% figure;imagesc(real(reshape(ffb,[Nx,Ny])));axis equal tight;  colormap(jet(2^13));
%lazy bydlocode
Nharmonics=50;
[mie_b,mie_c]=mie_coefs(k0,EF,eps,Nharmonics);

field_mie=(mie_cyl_die(R,EF,eps,Th+pi/2+alpha,mie_b,mie_c,k0));%-double(R<EF).*field_0;
%  +real(reshape(full(b),[2^dx,2^dy]))-double(R<EF).*field_0  double(R<EF).*

% figure('name','Mie series for PEC');imagesc(real(field_0));axis equal tight;  colormap(jet(2^13));colorbar;caxis([0,2.7])
figure('name','Mie series for PEC');imagesc(real(field_mie));axis equal tight;  colormap(jet(2^13));colorbar;