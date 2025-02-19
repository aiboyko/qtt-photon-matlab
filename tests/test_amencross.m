% test_amencross
%% 2D
R=1;

dr=15;
NR=2^dr;
hr=R/NR;

dphi=dr;
NPhi=2^dphi;
hphi=2*pi/NPhi;

tol=1e-5;%params.tol;

tt_idx_r=round(tt_x(2,dr)+1,1e-14);
tt_idx_phi=tt_x(2,dphi);

rcenter_handle = @(i) 2/3 *hr* ( (3*i.*(i - 1) + 1) ./ (2*i - 1) )*sin(hphi)/(hphi); %hardcoded fprmula for exact centers of masses of partial circle section that are our elements
area_i_handle = @(i) pi* hr^2 * (2*i-1)/NPhi;

tt_idx_rr=mtkron(tt_idx_r,tt_ones(2,dphi)); %correct!
tt_idx_phiphi=mtkron(tt_ones(2,dr),tt_idx_phi);

ij2rphi=@(ij) [rcenter_handle(ij(:,1))  ,  hphi*ij(:,2)];
rphi2xy=@(rphi)[rphi(:,1).*cos(rphi(:,2)) , rphi(:,1).*sin(rphi(:,2))];
planewave_h=@(xy) exp(1j*( params.rhs.k0x*xy(:,1)+params.rhs.k0y*xy(:,2) ));

rh_rfrf= @(xy)gh_xy_Nyst_selfheal(xy,params);

totalfun_gh= @(ij)rh_rfrf(rphi2xy(ij2rphi(ij)));
totalfun_gh3=@(ij)rh_rfrf(rphi2xy(ij2rphi(ij))).*area_i_handle(ij(:,1));

% piss=round(amen_cross({tt_idx_rr,tt_idx_phiphi},totalfun_gh3,tol*1e-1,'trunc_method','svd','kickrank',100, 'nswp',100,'max_err_jumps',9,'zrank',13),tol)
% figure; imagesc(real(reshape(full(piss),[NR,NPhi]))); colormap(redblue(2000));colorbar;


%% 4Darea_i_handle

% % % generating standard 4D mesh
idx_rrrr1=mtkron(tt_ones(2,dphi),tt_ones(2,dr), tt_ones(2,dphi), tt_idx_r); %correct!
idx_ffff1=mtkron(tt_ones(2,dphi),tt_ones(2,dr), tt_idx_phi,      tt_ones(2,dr));
idx_rrrr2=mtkron(tt_ones(2,dphi),tt_idx_r,      tt_ones(2,dphi), tt_ones(2,dr));
idx_pppp2=mtkron(tt_idx_phi     ,tt_ones(2,dr), tt_ones(2,dphi), tt_ones(2,dr));

% % % doing magical permute
idx_rrrr1=tt_standardpermute(idx_rrrr1,tol);
idx_ffff1=tt_standardpermute(idx_ffff1,tol);
idx_rrrr2=tt_standardpermute(idx_rrrr2,tol);
idx_pppp2=tt_standardpermute(idx_pppp2,tol);

ijij2rfrf=@(ijij)[rcenter_handle(ijij(:,1)), hphi*ijij(:,2), rcenter_handle(ijij(:,3)), hphi*ijij(:,4)];
ijij2frfr=@(ijij)[hphi*ijij(:,2),rcenter_handle(ijij(:,1)),  hphi*ijij(:,4),rcenter_handle(ijij(:,3))];
rfrf2xxyy=@(rf)[rf(:,1).*cos(rf(:,2)),rf(:,3).*cos(rf(:,4)) , rf(:,1).*sin(rf(:,2)), rf(:,3).*sin(rf(:,4))];
frfr2xxyy=@(rf)[rf(:,2).*cos(rf(:,1)),rf(:,4).*cos(rf(:,3)) , rf(:,2).*sin(rf(:,1)), rf(:,4).*sin(rf(:,3))];

% gh=@(xxyy) gh_xxyy_Nyst_selfheal(xxyy,params);
rh_rfrf=@(rfrf) sqrt(  rfrf(:,1).^2+rfrf(:,3).^2-2*rfrf(:,1).*rfrf(:,3).*cos( rfrf(:,2) - rfrf(:,4) )  );
pttv1=round(amen_cross({idx_rrrr1,idx_ffff1,idx_rrrr2,idx_pppp2},@(ijij)  (rh_rfrf(ijij2rfrf(ijij))),...
    tol*1e-1,'trunc_method','svd','kickrank',100, 'nswp',50,'max_err_jumps',5,'zrank',10),tol)


% gh=@(r) besselh( 0,2,2*pi*r/R );
gh=@(r) 1./(r+1e-10);
pttv2=round(amen_cross({pttv1},gh,...
    tol,'trunc_method','svd','kickrank',120, 'nswp',50,'max_err_jumps',20,'zrank',20),tol)

 G1=round(tt_vec2mat(pttv2),tol)
G2=round(tt_vec2mat(pttv2),tol);
% figure;imagesc(real( reshape(full(G1),[2^(dr+dphi),2^(dr+dphi)] ))); 
% colormap(redblue(1000)); 
% caxis([0,2]*R)
% colorbar
% axis equal tight



tic
% pttv2=round(amen_cross({idx_rrrr1,idx_ffff1,idx_rrrr2,idx_pppp2},@(ijij)  gh(frfr2xxyy(ijij2frfr(ijij))),...
%     tol*1e-1,'trunc_method','svd','kickrank',100, 'nswp',50,'max_err_jumps',5,'zrank',10),tol);

toc
