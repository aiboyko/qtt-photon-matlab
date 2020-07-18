clear all
close all
str=open('C:\Users\Alexey\Documents\qtt-photon\test_newfeatures\H from E\sergey_magneticfield.mat');


H=permute(str.H,[1 3 4 2]);
E=permute(str.E,[1 3 4 2]);
eps=permute(str.EPS,[2 1]);
[XX,YY]=meshgrid(str.Z,str.X);
Nx=max(size(str.Z));
Ny=max(size(str.X));
figure


subplot(1,2,1)
imagesc(XX(:,1),YY(1,:),real(E(:,:,1)));axis equal tight; colormap(redblue(1024)) ;colorbar ; caxis(2*[-1,1]);
subplot(1,2,2)
imagesc(XX(:,1),YY(1,:),real(E(:,:,3)));axis equal tight; colormap(redblue(1024));colorbar ; caxis(2*[-1,1]);


%Diff_test 
lambda=800;
hh=XX(1,2)-XX(1,1);

Dx=getdiffoper(Nx,1);
Dy=getdiffoper(Ny,1);
DDx=kron(Dx,eye(Ny)); %derivative to the right
DDy=kron(eye(Nx),Dy); % derivative downwards

% vx=([1:Nx]'.^2)/2;
% vy=[1:Ny]';
% vv=kron(vx,vy); %vv may as well be a picture

vv=H(:,:,2);
vv=vv(:);

dvvx=(-1j*lambda/(2*pi*hh)*DDx*vv);
dvvx=reshape(dvvx,[Ny,Nx])./eps;

dvvy=( 1j*lambda/(2*pi*hh)*DDy*vv);
dvvy=reshape(dvvy,[Ny,Nx])./eps;

figure;
subplot(1,2,1)
imagesc(real(dvvx));axis equal tight;colormap(redblue(1024));colorbar;caxis(2*[-1,1])
subplot(1,2,2)
imagesc(real(dvvy));axis equal tight;colormap(redblue(1024));colorbar;caxis(2*[-1,1])
