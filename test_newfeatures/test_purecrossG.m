%Purecross implemented successfully!

% close all
% dx=4;
% dy=4;
% d=dx+dy;
% Nx=2^dx;
% Ny=2^dy;
% tol=1e-5;%1e-1*2^-dx
% softening_coef=1e-3;
% Lx=1;
% Ly=1;%2*pi;
% 
% k0=2*pi/Lx;
% 
% hx=Lx/(Nx);
% hy=Lx/(Nx);
% 
ttX=hx/2+hx*tt_x(2,dx);%volume-centric lattice is assumed
ttY=hy/2+hy*tt_x(2,dy);
% 
% params.hx=hx;
% params.hy=hy;
% params.dx=dx;
% params.dy=dy;
% params.tol=tol;
% params.k0=k0;

%general approach with composition of functions
 %gh=@(xxyy) exp(1j*k0*sqrt(   (xxyy(:,1)-xxyy(:,2)).^2+(xxyy(:,3)-xxyy(:,4)).^2)  )./sqrt( (xxyy(:,1)-xxyy(:,2)).^2+(xxyy(:,3)-xxyy(:,4)).^2+tol*1e+1)
 
gh=@(xxyy)-hx*hy* 1j/4* besselh( 0,1,k0*hx*softeningcoef+k0*sqrt(  (xxyy(:,1)-xxyy(:,2)).^2+(xxyy(:,3)-xxyy(:,4)).^2 ));
gh_curv2xy=@(coord) [coord(:,1).*cos(coord(:,3)),coord(:,2).*cos(coord(:,4)), coord(:,1).*sin(coord(:,3)) , coord(:,2).*sin(coord(:,4))]

%direct computation woth cosine theorem
%gh_total=@(coord) 1./sqrt( coord(:,1).^2+coord(:,2).^2 - 2*coord(:,1).*coord(:,2).*cos(coord(:,3)-coord(:,4))  +   tol*1e1)

ttX1=mtkron(ttX,tt_ones(2,dy),tt_ones(2,dx),tt_ones(2,dy)); %correct!
ttY1=mtkron(tt_ones(2,dx),ttY,tt_ones(2,dx),tt_ones(2,dy));
ttX2=mtkron(tt_ones(2,dx),tt_ones(2,dy),ttX,tt_ones(2,dy));
ttY2=mtkron(tt_ones(2,dx),tt_ones(2,dx),tt_ones(2,dy),ttY);

pttX1=tt_standardpermute(ttX1,tol);
pttX2=tt_standardpermute(ttX2,tol);
pttY1=tt_standardpermute(ttY1,tol);
pttY2=tt_standardpermute(ttY2,tol);

% pttv=round(amen_cross({pttX1,pttX2,pttY1,pttY2},@(x)gh(gh_curv2xy(x)),tol,...
% 'trunc_method','svd','kickrank',100, 'nswp',50,'max_err_jumps',3,'zrank',10),tol);

pttv=round(amen_cross({pttX1,pttX2,pttY1,pttY2},@(x)gh(x),tol,...
'trunc_method','svd','kickrank',100, 'nswp',50,'max_err_jumps',3,'zrank',10),tol);
modes = pttv.n(1:d);
ttm_custom=round(tt_vec2mat(pttv),tol)
% maxmrank=max(rank(ttm_custom))
figure('name','ttmatrix with custom reshape');imagesc(imag(full(ttm_custom)));colormap(jet(2^12));colorbar;axis equal tight

patch_for_d4=0;
gh2=@(xxyy) gh_xxyy_LCN(xxyy,patch_for_d4,params);
pttv2=round(amen_cross({pttX1,pttX2,pttY1,pttY2},@(x)gh2(x),tol,...
'trunc_method','svd','kickrank',100, 'nswp',50,'max_err_jumps',3,'zrank',10),tol);
modes = pttv2.n(1:d);
ttm_custom2=round(tt_vec2mat(pttv2),tol)
figure('name','ttmatrix with custom reshape and LCN function handle');imagesc(abs(full(ttm_custom2)));colormap(jet(2^12));colorbar;axis equal tight



% maxvrank=max(rank(pttv))

%full
% [X,Y]=meshgrid(0:Nx-1,0:Ny-1);
% X1=full(ttX1);
% X2=full(ttX2);
% Y1=full(ttY1);
% Y2=full(ttY2);
% honest_v=gh([X1,X2,Y1,Y2]);
% honest_ttv=tt_tensor(honest_v,tol,2*ones(1,2*d))
% honest_m=reshape(full(honest_ttv),2^d*[1 1]);
% figure('name','honestly calculated matrix');imagesc(honest_m);colormap(jet(2^12));colorbar;axis equal tight
% caxis([0,20])



% corrected_diag=666;
% corrected_ttm=ttm_custom+corrected_diag*tt_eye(2,d)-diag(diag(ttm_custom));

