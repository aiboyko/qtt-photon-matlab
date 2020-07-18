%test1
R=2;

% creating RPHI-mesh

Rmax = 1;
NR=50;
Nphi=NR;
R=Rmax*(0:Nphi-1)/(Nphi-1);
PHI=2*pi*(0:Nphi-1)/Nphi;
[RR,PHIPHI]=meshgrid(R,PHI);

%creating XY-mesh
xleft=-Rmax;
xright=Rmax;
yleft=-Rmax;
yright=Rmax;
Nx=NR;
Ny=NR;
Lx=xright-xleft;
Ly=yright-yleft;
hx=Lx/Nx;
hy=Ly/Ny;
X=xleft + hx/2 :hx:xright-hx/2;
Y=yleft + hy/2 :hy:yright-hy/2;

[XXq,YYq]=meshgrid(X,Y);

%creating data in RPHI coordinates

polardata= zeros(NR);
polardata(1:ceil(NR/3),1:ceil(NR/2))=ones(ceil(NR/3),ceil(NR/2));
% 
% [XX,YY] = pol2cart(PHIPHI,RR);
% 
% 
% tic
% im_out=griddata(XX,YY,polardata,XXq,YYq,'natural');axis equal tight
% toc

tic
RRq=sqrt(XXq.^2+YYq.^2);
PHIPHIq = angle(XXq+1j*YYq);

poldata_out=interp2(RR,PHIPHI,polardata,RRq,PHIPHIq);
toc

figure('name','interp2');
imagesc(poldata_out);
axis equal tight;
colorbar


figure('name','scatinterp');
imagesc(im_out);
axis equal tight;
colorbar