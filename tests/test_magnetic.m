%test_ pretend we were solving helmholtz for B, then how to recover Ex, Ez:



%Diff_test 

i=1
Dx=getdiffoper(params.grid.Nx,1);
Dy=getdiffoper(params.grid.Ny,1);
DDx=kron(Dx,eye(params.grid.Ny)); %derivative to the right
DDy=kron(eye(params.grid.Nx),Dy); % derivative downwards

% vx=([1:Nx]'.^2)/2;
% vy=[1:Ny]';
% vv=kron(vx,vy); %vv may as well be a picture

 %soo we obtain Ex,Ez from Hy
vv=data_fft{1};
vv=vv(:);
% vv=permute(reshape(full(x_qtt),[params.grid.Nx,params.grid.Ny]),[2 1]);
% vv=vv(:);

% ffmask=full(material_mask);

dvvx=(-1j*params.rhs.lambda/(2*pi*params.grid.hx)*DDx*vv);
dvvx=reshape(dvvx./(1+ffmask(:)*(params.parts.eps-1)),[params.grid.Ny,params.grid.Nx]);

dvvy=( 1j*params.rhs.lambda/(2*pi*params.grid.hy)*DDy*vv);
dvvy=reshape(dvvy./(1+ffmask(:)*(params.parts.eps-1)),[params.grid.Ny,params.grid.Nx]);

figure;
subplot(1,2,1)
imagesc(imag(dvvx));axis equal tight;colormap(redblue(1024));colorbar;
 caxis(2*[-1,1])
subplot(1,2,2)
imagesc(imag(dvvy));axis equal tight;colormap(redblue(1024));colorbar;
 caxis(2*[-1,1])
