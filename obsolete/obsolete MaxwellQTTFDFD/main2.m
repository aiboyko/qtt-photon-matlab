%By Alexey Boyko, Skolkovo Institute of Science and Technology,
%FDFD QTT Photon solver
%a.boyko@skoltech.ru
%alexey.boyko@skolkovotech.ru
%Moscow, Russia, 2015


%weirdly enough, only if one choses the tf\sf boundary at PML boundary,
%everything works

clear all;
%% Don't touch that!
CONST.AXIS.X=1;
CONST.AXIS.Y=2;
CONST.AXIS.Z=3;
CONST.FIELD.E=1;
CONST.FIELD.H=2;
%% Setting up general constants -------wtha------------------------
tol=1e-7;
SparseSolverON=1;
QTTsolverON=0;
inc_angle=0/360*2*pi;
%% Setting up the grid parameters------------------------------

dx=7;dy=7;dz=1;
d=dx+dy+dz;

doUPML=0;

tfsffrac=0.30;

XisPeriodic=0;
YisPeriodic=1;
ZisPeriodic=0;

%point for probing magnitude of reflected wave
 xtestidx=floor(2^dx*tfsffrac)-2;
 ytestidx=2^(dy-1);
 ztestidx=2^(dz-1);

%% Setting up physical envinronment 
eps=1;

CONST_UNITS={'SI','Gauss'};

units='Gauss'
params.units=units;

switch units
case 'SI'
    Lx=5e-6;
    Ly=Lx;
    Lz=Lx/(2^dx-1)*(2^dz-1);
    mu=1.0;
    eps0=8.85e-12;
    mu0=4*pi*1e-7;
    C=3e8;
    lambda=500e-9;%2^(dy-4)+2^(dy-8);
    k0=2*pi/lambda;
    omega=k0*C;
    etha0=sqrt(mu0/eps0);
case 'Gauss'       
    Lx=5e-4;
    Ly=Lx;
    Lz=Lx/(2^dx-1)*(2^dz-1);
    mu=1.0;
    eps0=1;
    mu0=1;
    C=3e10;
    lambda=500e-7;%2^(dy-4)+2^(dy-8);
    k0=2*pi/lambda;
    omega=k0*C;
    etha0=1;
end

hx=Lx/(2^dx-1); 
hy=Ly/(2^dy-1);
hz=Lz/(2^dz-1);

%% renormalization of length. should work for gauss+wave, dunno about gauss+EH, shouldn't work for SI
doRenorm=1
if doRenorm
    hmin=min([hx,hy,hz]);
    Lx=Lx/hmin;
    Ly=Ly/hmin;
    Lz=Lz/hmin;

    %our unit of length is now this
    k0=k0*hmin; %renormalized
    hx=hx/hmin;
    hy=hy/hmin;
    hz=hz/hmin;
    C=C/hmin;
    lambda=lambda/hmin;

    if strcmp(units,'SI')
        eps0=eps0*hmin;
        mu0=mu0*hmin;
    end
end
%% 

xtest=xtestidx*hx;
ytest=ytestidx*hy;
ztest=ztestidx*hz;

k0x=k0*cos(inc_angle);
k0y=k0*sin(inc_angle);

BlochBCS_EX=exp(1j*k0x*(Lx));
BlochBCS_HX=BlochBCS_EX';
BlochBCS_EY=exp(1j*k0y*(Ly));
BlochBCS_HY=BlochBCS_EY';


%% Setting up geometry and calculating in on qtt lattice-------
hole=2*lambda;
holes=2^(dy-4);
depth=holes;

params.k0=k0;

params.tol=tol;
params.eps0=eps0;
params.mu0=mu0;

params.dx=dx;
params.dy=dy;
params.dz=dz;

params.Lx=Lx;
params.Ly=Ly;
params.Lz=Lz;

params.hx=hx;
params.hy=hy;
params.hz=hz;

params.omega=omega;
params.etha0=etha0;
pml.thickness=0.8*lambda;

pml.xlowwidth=pml.thickness;
pml.xhighwidth=pml.thickness;
pml.ylowwidth=pml.thickness;
pml.yhighwidth=pml.thickness;
pml.zlowwidth=0;
pml.zhighwidth=0;
pml.R=1e-7;
pml.m=3;

if doUPML == 1
    [syz_x,sxz_y,sxy_z,invsxz_y,invsyz_x,invsxy_z]=createUPML(params,pml);
    
else
      [invsx,invsy,invsz]=createSCPML2(params, pml);
    if SparseSolverON
        spdiag_invsx=diag(sparse(full(invsx)));
        spdiag_invsy=diag(sparse(full(invsy)));
        spdiag_invsz=diag(sparse(full(invsz)));
    end
    if QTTsolverON
        diag_invsx=diag(invsx);
        diag_invsy=diag(invsy);
        diag_invsz=diag(invsz);
    end
end

%both of index boundaries are included(!) in the box object
% xb1_1=2^(dx-1)-floor(depth/2);xb1_2=2^(dx-1)+floor(depth/2);
% yb1_1=1;yb1_2=2^(dy);
% zb1_1=1;zb1_2=2^dz;
% box1=geom_box(dx,dy,dz,tol,xb1_1,xb1_2,yb1_1,yb1_2,zb1_1,zb1_2);

geom={};

xb1=2^(dx-1)-floor(depth/2);
xb2=2^(dx-1)+floor(depth/2);

zb1=1;
zb2=2^dz;

mask_of_coords=tt_zeros([2*ones(1,d)]);
nboxmax=(2^dy)/holes;
% for i=1:2:nboxmax
% geom{i}=geom_box(params,xb1,xb2,i*holes-floor(holes/2),i*holes+floor(holes/2)-1,zb1,zb2);
% mask_of_coords=round(mask_of_coords+geom{i},tol);
% end
geom{nboxmax+1}=geom_box(params,xb1-depth,xb1-1,1,2^dy,zb1,zb2);

mask_of_coords=round(mask_of_coords+geom{nboxmax+1},tol);

% general=@(xyz) exp(2*pi*1j*xyz(1)/lambda);
% mask_amen=amen_cross(2*ones(d, 1),@(bin) general(bin2xyz(bin-1,...
%      params)), 1e-14,'nswp',40)
if d<=23 && SparseSolverON
[fullmask1_3D, X, Y, Z]  =  F1_LtoF1_3D(full(mask_of_coords), params);
end
% figure();imagesc(fullmask1_3D(:,:,2));colorbar
% params.mask3d=fullmask1_3D;
% 
letmePlot3d=0;
if letmePlot3d
figure();
iso_mask=isosurface(X,Y,Z,fullmask1_3D);
p=patch(iso_mask);
p.FaceAlpha=0.9;
p.EdgeColor   =   'none';
daspect([1,1,1])
view(3); axis tight
p.FaceColor   =   'yellow';
camlight 
lighting gouraud
end

% 
% xb1_1=2^(dx-1)-floor(depth/2);xb1_2=2^(dx-1)+floor(depth/2);
% yb1_1=1;yb1_2=2^(dy);
% zb1_1=1;zb1_2=2^dz;
% box=geom_box(dx,dy,dz,tol,xb1_1,xb1_2,yb1_1,yb1_2,zb1_1,zb1_2);
% 

%scattered box
xb3_1=floor(2^dx*tfsffrac);%ceil(pml.xlowwidth/hx)+1 ;
xb3_2=2^dx%-ceil(pml.xlowwidth/hx);
yb3_1=1;%ceil(pml.ylowwidth/hy)+1 ;
yb3_2=2^dy;%2^dy-ceil(pml.ylowwidth/hy);
zb3_1=1;
zb3_2=2^dz;

box_Q=round(geom_box(params,xb3_1,xb3_2,yb3_1,yb3_2,zb3_1,zb3_2),tol);
% box_Q_m=diag(box_Q);
% nil=ttm_zeros(2,dx+dy+dz);
% box_Q_All=round(tt_blockwise({box_Q_m,      nil,       nil ;
%                             nil,           box_Q_m,  nil;
%                             nil,           nil,       box_Q_m     }),tol)

if SparseSolverON 
spbox_Q=sparse(full(box_Q));
 Q=diag([spbox_Q;spbox_Q;spbox_Q]);
 Q_EH=diag([spbox_Q;spbox_Q;spbox_Q;spbox_Q;spbox_Q;spbox_Q]);
end
if QTTsolverON
 box_3Q=round(tt_blockwise({box_Q;box_Q;box_Q}),tol);
 box_6Q_EH=round(tt_blockwise({box_3Q;box_3Q}),tol);
 ttQ=diag(box_3Q);
 ttQ_EH=diag(box_6Q_EH);
%  spbox_Q=speye(size(spbox_Q))-spbox_Q
end

%tt_ones(2,d)% round(box1+box2,tol);
%% source creation
N=2^d;
% planewave1_3D  =  @(x)exp(-1j*k0*x ).* double(x<200);
% planewave1_3D  =  @(x)exp(-1j*k0*x ).* double(x>2*pml_xlow&(x<(2^dx-2*pml_xhigh-1)));
 planewave1_3D  =  @(x,y) exp(1j*(k0x*x+k0y*y));
 if (d<=22)&&SparseSolverON
     SOURCEFIELD1_3D=planewave1_3D(X,Y);
     SOURCEFIELD3_L=[zeros(N,1);zeros(N,1);F1_3DtoF1_L(SOURCEFIELD1_3D)];
     SOURCEFIELD3_L_EH=[SOURCEFIELD3_L;zeros(N,1);1/etha0*F1_3DtoF1_L(SOURCEFIELD1_3D);zeros(N,1)];
 end

%% Assemblying Maxwell operator full\sparse format------------------------------
if SparseSolverON
nil_f=sparse(2^d,2^d);

mask_of_coords_f=sparse(round(full(mask_of_coords)));

%here there should be a generalization for complex geometry
eps_3_L=mask_of_coords_f*eps + (1-mask_of_coords_f);

if doUPML
    sp_x=sparse(full(syz_x));
    sp_y=sparse(full(sxz_y));
    sp_z=sparse(full(sxy_z));

    eps_Ex_f=sparse(eps_3_L.*sp_x);
    eps_Ey_f=sparse(eps_3_L.*sp_y);
    eps_Ez_f=sparse(eps_3_L.*sp_z);

    eps_Ex_m_f=sparse(diag(eps_Ex_f));
    eps_Ey_m_f=sparse(diag(eps_Ey_f));
    eps_Ez_m_f=sparse(diag(eps_Ez_f));

    eps_All_f=[      eps_Ex_m_f,    nil_f,        nil_f;
                     nil_f,         eps_Ey_m_f,   nil_f;
                     nil_f,         nil_f,        eps_Ez_m_f     ];

    invspyz_x=sparse(full(invsyz_x));
    invspxz_y=sparse(full(invsxz_y));
    invspxy_z=sparse(full(invsxy_z));

    inv_eps_3_L=mask_of_coords_f*eps^(-1) + (1-mask_of_coords_f);

    inv_eps_Ex_f=sparse(inv_eps_3_L.*invspyz_x);
    inv_eps_Ey_f=sparse(inv_eps_3_L.*invspxz_y);
    inv_eps_Ez_f=sparse(inv_eps_3_L.*invspxy_z);

    inv_eps_Ex_m_f=diag(inv_eps_Ex_f);
    inv_eps_Ey_m_f=diag(inv_eps_Ey_f);
    inv_eps_Ez_m_f=diag(inv_eps_Ez_f);

    inv_eps_All_f=[inv_eps_Ex_m_f, nil_f,          nil_f;
                nil_f,         inv_eps_Ey_m_f,  nil_f;
                nil_f,         nil_f,          inv_eps_Ez_m_f     ];

    %mu goes here
    %mu_3_L === 1
    mu_Ex_f=sparse(sp_x);
    mu_Ey_f=sparse(sp_y);
    mu_Ez_f=sparse(sp_z);

    mu_Ex_m_f=sparse(diag(mu_Ex_f));
    mu_Ey_m_f=sparse(diag(mu_Ey_f));
    mu_Ez_m_f=sparse(diag(mu_Ez_f));

    mu_All_f=[      mu_Ex_m_f,    nil_f,       nil_f;
                    nil_f,        mu_Ey_m_f,   nil_f;
                    nil_f,        nil_f,       mu_Ez_m_f     ];
            
    inv_mu_Ex_f=sparse(mask_of_coords_f*mu^(-1) + (1-mask_of_coords_f));
    inv_mu_Ey_f=sparse(mask_of_coords_f*mu^(-1) + (1-mask_of_coords_f));
    inv_mu_Ez_f=sparse(mask_of_coords_f*mu^(-1) + (1-mask_of_coords_f));

    inv_mu_Ex_m_f=diag(inv_mu_Ex_f.*invspyz_x);
    inv_mu_Ey_m_f=diag(inv_mu_Ey_f.*invspxz_y);
    inv_mu_Ez_m_f=diag(inv_mu_Ez_f.*invspxy_z);


    inv_mu_All_f=[inv_mu_Ex_m_f, nil_f,        nil_f;
                nil_f,         inv_mu_Ey_m_f,  nil_f;
                nil_f,         nil_f,          inv_mu_Ez_m_f     ];
else
    eps_Ex_f=sparse(eps_3_L);
    eps_Ey_f=sparse(eps_3_L);
    eps_Ez_f=sparse(eps_3_L);

    eps_Ex_m_f=sparse(diag(eps_Ex_f));
    eps_Ey_m_f=sparse(diag(eps_Ey_f));
    eps_Ez_m_f=sparse(diag(eps_Ez_f));

    eps_All_f=[      eps_Ex_m_f,    nil_f,        nil_f;
                     nil_f,         eps_Ey_m_f,   nil_f;
                     nil_f,         nil_f,        eps_Ez_m_f     ];

    inv_mu_Ex_f=sparse(mask_of_coords_f*mu^(-1) + (1-mask_of_coords_f));
    inv_mu_Ey_f=sparse(mask_of_coords_f*mu^(-1) + (1-mask_of_coords_f));
    inv_mu_Ez_f=sparse(mask_of_coords_f*mu^(-1) + (1-mask_of_coords_f));

    inv_mu_Ex_m_f=diag(inv_mu_Ex_f);
    inv_mu_Ey_m_f=diag(inv_mu_Ey_f);
    inv_mu_Ez_m_f=diag(inv_mu_Ez_f);

    inv_mu_All_f=[inv_mu_Ex_m_f, nil_f,        nil_f;
                nil_f,         inv_mu_Ey_m_f,  nil_f;
                nil_f,         nil_f,          inv_mu_Ez_m_f     ];             

end
%warning! mtkron takes tensor product is reverse order
%for vector wave equation 

DEX=sparse(diag(ones(2^dx-1,1),1)-diag(ones(2^dx,1),0));
DEX(2^dx,1)=BlochBCS_EX*XisPeriodic; %periodic 
DEY=sparse(diag(ones(2^dy-1,1),1)-diag(ones(2^dy,1),0));
DEY(2^dy,1)=BlochBCS_EY*YisPeriodic; %periodic 
DEZ=sparse(diag(ones(2^dz-1,1),1)-diag(ones(2^dz,1),0));
DEZ(2^dz,1)=ZisPeriodic; %periodic 
%so right now X is the largest harmonic, Z is the highest-freq harmonic
DHX=sparse(-diag(ones(2^dx-1,1),-1)+diag(ones(2^dx,1),0));
DHX(1,2^dx) = -XisPeriodic*BlochBCS_HX;
DHY=sparse(-diag(ones(2^dy-1,1),-1)+diag(ones(2^dy,1),0));
DHY(1,2^dy) = -YisPeriodic*BlochBCS_HY;
DHZ=sparse(-diag(ones(2^dz-1,1),-1)+diag(ones(2^dz,1),0));
DHZ(1,2^dz) = -ZisPeriodic;

deltaEx_f=kron(1/hx * DEX,   kron( speye(2^dy) ,     speye(2^dz)));
deltaEy_f=kron(speye(2^dx),  kron( 1/hy * DEY,       speye(2^dz)));
deltaEz_f=kron(speye(2^dx),  kron( speye(2^dy) ,     1/hz * DEZ ));
deltaHx_f=kron(1/hx * DHX, kron( speye(2^dy) ,     speye(2^dz)  ) );
deltaHy_f=kron(speye(2^dx),  kron( 1/hy * DHY,     speye(2^dz)));
deltaHz_f=kron(speye(2^dx),  kron( speye(2^dy) ,     1/hz * DHZ));

curlE_f=[ nil_f,        -deltaEz_f,       deltaEy_f ;
          deltaEz_f,        nil_f,            -deltaEx_f;
         -deltaEy_f,       deltaEx_f,        nil_f     ];

curlH_f=[nil_f,           -deltaHz_f,       deltaHy_f ;
         deltaHz_f,        nil_f,            -deltaHx_f;
         -deltaHy_f,       deltaHx_f,        nil_f     ]; 
     
if ~doUPML
    sc_deltaEx_f=spdiag_invsx*deltaEx_f;
    sc_deltaHx_f=spdiag_invsx*deltaHx_f;
    sc_deltaEy_f=spdiag_invsy*deltaEy_f;
    sc_deltaHy_f=spdiag_invsy*deltaHy_f;
    sc_deltaEz_f=spdiag_invsz*deltaEz_f;
    sc_deltaHz_f=spdiag_invsz*deltaHz_f;
    
    
    sc_curlE_f=[ nil_f,        -sc_deltaEz_f,       sc_deltaEy_f ;
          sc_deltaEz_f,        nil_f,            -sc_deltaEx_f;
         -sc_deltaEy_f,       sc_deltaEx_f,        nil_f     ];

    sc_curlH_f=[nil_f,           -sc_deltaHz_f,       sc_deltaHy_f ;
         sc_deltaHz_f,        nil_f,            -sc_deltaHx_f;
         -sc_deltaHy_f,       sc_deltaHx_f,        nil_f     ]; 
end
% compressing deltas doesnt help much


% coefs_empty_pml=diag([sp_x;sp_y;sp_z]);
% coefs_inv_empty_pml=diag([sparse(invspyz_x);sparse(invspxz_y);sparse(invspxy_z)]);

if doUPML
    switch(units)
    case 'SI'
        Aff=curlH_f*inv_mu_All_f*curlE_f/mu0-omega^2*eps_All_f*eps0;
        A_EH_ff=[1j*omega*eps0*eps_All_f,-curlH_f;...
        curlE_f, 1j*omega*mu0*mu_All_f];
    case 'Gauss'
        Aff=curlH_f*inv_mu_All_f*curlE_f-k0^2*eps_All_f;
        A_EH_ff=[1j*k0*eps_All_f,-curlH_f;...
        sc_curlE_f, 1j*k0*mu0*speye(size(eps_All_f))];
    end
else
    switch(units)
    case 'Gauss'
        Aff=sc_curlH_f*inv_mu_All_f*sc_curlE_f-k0^2*eps_All_f;
        A_EH_ff=[1j*k0*eps_All_f,-sc_curlH_f;...
        sc_curlE_f, 1j*k0*mu0*speye(size(eps_All_f))];
    case 'SI'
        Aff=sc_curlH_f*inv_mu_All_f*sc_curlE_f/mu0-omega^2*eps_All_f*eps0;
        A_EH_ff=[1j*omega*eps0*eps_All_f,-sc_curlH_f;...
        sc_curlE_f, 1j*omega*mu0*speye(size(eps_All_f))];
    end
end
% A_empty_ff=curlH_f*curlE_f-k0^2*speye(2^d * 3) ;
% A_empty_ff_scpml=sc_curlH_f*sc_curlE_f-k0^2*speye(2^d * 3) ;


% 
% A_empty_EH_ff=[1j*k0*speye(size(eps_All_f)), -sc_curlH_f;...
%                 sc_curlE_f, 1j*k0*speye(size(eps_All_f))];

% A_empty_ff_upml=curlH_f*coefs_inv_empty_pml*curlE_f-k0^2*coefs_empty_pml;

% P=inv_eps_All_f;

b_ff=(Q*Aff-Aff*Q)*SOURCEFIELD3_L;
b_EH_ff=(Q_EH*A_EH_ff-A_EH_ff*Q_EH)*SOURCEFIELD3_L_EH;

SparseKrylovON=0;
if SparseSolverON && SparseKrylovON
    sparseiter=23;
[F3_L_krylov, flag, relres,iter,reshist]=gmres(Aff, b_ff,sparseiter,...
    1e-14,sparseiter,[],[],SOURCEFIELD3_L);
figure('Name','Vector wave Krylov (bicgstab) convergence');
plot(log10(reshist));
[F3_L_EH_krylov, flag2, relres2,iter2,reshist2]=gmres(A_EH_ff, b_EH_ff,sparseiter,...
    1e-14,sparseiter,[],[],SOURCEFIELD3_L_EH);
figure('Name','EH Krylov (bicgstab) convergence');
plot(log10(reshist2));

Fz_s_3d_krylov=F1_LtoF1_3D(F3_L_krylov((2*N+1):(3*N)),params);
Fz_s_3d_EH_krylov=F1_LtoF1_3D(F3_L_EH_krylov((2*N+1):(3*N)),params);

figure('Name','Krylov - full - Vector Wave');imagesc(real(Fz_s_3d_krylov(:,:,2^(dz-1))));
title('Ez')
axis equal tight
colorbar;

figure('Name','Krylov - full - EH');imagesc(real(Fz_s_3d_EH_krylov(:,:,2^(dz-1))));
title('Ez')
axis equal tight
colorbar;

end
% flag
% bicgstab_res=norm(Aff*F3_L-bf3_L,1)

% Direct solver for sparse
% F3_L_source_pml=A_empty_ff_pml\bf3_L;
% bf3_L_1=A_empty_ff_pml
% F3_L_orig=Aff\bf3_L;

F3_L=Aff\b_ff;
F3_L_EH=A_EH_ff\b_EH_ff;

ftt=tt_tensor( F3_L,tol,[2*ones(1,d),3]);
F3_L=full(ftt);
F3_L_3d =    F1_LtoF1_3D(F3_L((2*N+1):(3*N)),params);
Fz_s_3d_EH = F1_LtoF1_3D(F3_L_EH((2*N+1):(3*N)),params);

F3L_visualize(F3_L_EH,params,'field','Ez','fig',1,'name','Direct solver - EH','func','re');
F3L_visualize(F3_L_EH,params,'field','Ez','fig',2,'name','Direct solver - EH','func','im');
% F3L_visualize(F3_L_EH,params,'field','Ey','fig',3,'name','Direct solver - EH');
% F3L_visualize(F3_L_EH,params,'field','Hx','fig',4,'name','Direct solver - EH');
% F3L_visualize(F3_L_EH,params,'field','Hy','fig',5,'name','Direct solver - EH');
% F3L_visualize(F3_L_EH,params,'field','Hz','fig',6,'name','Direct solver - EH');

EH_refl=abs(Fz_s_3d_EH(xtestidx+1,ytestidx+1,ztestidx+1))
% 
% F3L_visualize(F3_L,params,'field','Ex','fig',10,'name','Direct solver - Vector Wave');
% F3L_visualize(F3_L,params,'field','Ey','fig',11,'name','Direct solver - Vector Wave');
F3L_visualize(F3_L,params,'field','Ez','fig',12,'name','Direct solver - Vector Wave');

vectorwave_refl=abs(F3_L_3d(xtestidx+1,ytestidx+1,ztestidx+1))
% hold on
% contour(fullmask1_3D(:,:,2^(dz-1)),1,'color','Black');
% hold off


% dir_res=norm(Aff*F3_L1-bf3_L,1)
% figKr=plot3fields(F3_L,bf3_L,2^(dz-1), dx, dy, dz, hx, hy, hz,1,'Krylov Solver');

N=2^(dx+dy+dz);

% Fz_s_3d_orig = F1_LtoF1_3D(F3_L_orig((2*N+1):(3*N)),dx,dy,dz,hx,hy,hz);
% figure();imagesc(real(Fz_s_3d_orig(:,:,2^(dz-1))));
% title('Ez')
% axis equal tight
% colorbar;
% 
%  figDir=plot3fields(F3_L,bf3_L,fullmask1_3D, 2^(dz-1), params,2,'Direct Solver')

end
%% Assemblying Maxwell operator it qtt format------------------------------
if QTTsolverON
nil=ttm_zeros(2,dx+dy+dz);

eps_Ex=round((mask_of_coords*eps + (1-mask_of_coords)),tol);
eps_Ey=round((mask_of_coords*eps + (1-mask_of_coords)),tol);
eps_Ez=round((mask_of_coords*eps + (1-mask_of_coords)),tol);

eps_Ex_m=diag(eps_Ex);
eps_Ey_m=diag(eps_Ey);
eps_Ez_m=diag(eps_Ez);

eps_All=round(tt_blockwise({eps_Ex_m,      nil,       nil ;
                            nil,           eps_Ey_m,  nil;
                            nil,           nil,       eps_Ez_m     }),tol)

%--for preconditioner----------------------------
                        
% inv_eps_Ex=round(mask_of_coords*eps^(-1) + (1-mask_of_coords),tol);
% inv_eps_Ey=round(mask_of_coords*eps^(-1) + (1-mask_of_coords),tol);
% inv_eps_Ez=round(mask_of_coords*eps^(-1) + (1-mask_of_coords),tol);
% 
% inv_eps_Ex_m  =  diag(inv_eps_Ex);
% inv_eps_Ey_m  =  diag(inv_eps_Ey);
% inv_eps_Ez_m  =  diag(inv_eps_Ez);
% 
% inv_eps_All=round(tt_blockwise({inv_eps_Ex_m, nil,          nil;
%                                nil,         inv_eps_Ey_m,  nil;
%                                nil,         nil,          inv_eps_Ez_m     }),tol)           
% Ptt=inv_eps_All
%------------------------------------------------
                        
inv_mu_Ex=round((mask_of_coords*mu^(-1) + (1-mask_of_coords)),tol);
inv_mu_Ey=round((mask_of_coords*mu^(-1) + (1-mask_of_coords)),tol);
inv_mu_Ez=round((mask_of_coords*mu^(-1) + (1-mask_of_coords)),tol);

inv_mu_Ex_m  =  diag(inv_mu_Ex);
inv_mu_Ey_m  =  diag(inv_mu_Ey);
inv_mu_Ez_m  =  diag(inv_mu_Ez);

inv_mu_All=round(tt_blockwise({inv_mu_Ex_m, nil,          nil;
                               nil,         inv_mu_Ey_m,  nil;
                               nil,         nil,          inv_mu_Ez_m     }),tol);             

%warning! mtkron takes tensor product is reverse order
deltaEx=mtkron(tt_eye(2,dz),                tt_eye(2,dy),                   1/hx * fd_tt_DE(2,dx,XisPeriodic,tol,BlochBCS_EX));
deltaEy=mtkron(tt_eye(2,dz),                1/hy * fd_tt_DE(2,dy,YisPeriodic,tol,BlochBCS_EY),    tt_eye(2,dx));
deltaEz=mtkron(1/hz * fd_tt_DE(2,dz,ZisPeriodic,tol,1), tt_eye(2,dy),                   tt_eye(2,dx));

%so right now X is the largest harmonic, Z is the highest-freq harmonic

deltaHx=mtkron(tt_eye(2,dz),                tt_eye(2,dy),               1/hx * fd_tt_DH(2,dx,XisPeriodic,tol,BlochBCS_HX));
deltaHy=mtkron(tt_eye(2,dz),                1/hy * fd_tt_DH(2,dy,YisPeriodic,tol,BlochBCS_HY),tt_eye(2,dx));
deltaHz=mtkron(1/hz * fd_tt_DH(2,dz,ZisPeriodic,tol,1), tt_eye(2,dy),               tt_eye(2,dx));


sc_deltaEx=diag(invsx)*deltaEx;
sc_deltaEy=diag(invsy)*deltaEy;
sc_deltaEz=diag(invsz)*deltaEz;

sc_deltaHx=diag(invsx)*deltaHx;
sc_deltaHy=diag(invsy)*deltaHy;
sc_deltaHz=diag(invsz)*deltaHz;
% compressing deltas doesnt help much

% curlE=tt_blockwise({nil,           -deltaEz,       deltaEy ;
%                     deltaEz,        nil,            -deltaEx;
%                     -deltaEy,       deltaEx,        nil     });
% 
% curlH=tt_blockwise({nil,           -deltaHz,       deltaHy ;
%                     deltaHz,        nil,            -deltaHx;
%                     -deltaHy,       deltaHx,        nil     });   

                
sc_curlE=round(tt_blockwise({nil,           -sc_deltaEz,       sc_deltaEy ;
                    sc_deltaEz,        nil,            -sc_deltaEx;
                    -sc_deltaEy,       sc_deltaEx,        nil     }),tol);

sc_curlH=round(tt_blockwise({nil,           -sc_deltaHz,       sc_deltaHy ;
                    sc_deltaHz,        nil,            -sc_deltaHx;
                    -sc_deltaHy,       sc_deltaHx,        nil     }),tol);  
                          


uno=round(tt_blockwise({tt_eye(2,d),     nil,      nil;
                        nil,     tt_eye(2,d),      nil;
                        nil,     nil,      tt_eye(2,d)}),tol);
                    
nil3=round(tt_blockwise({nil,     nil,      nil;
                        nil,     nil,      nil;
                        nil,     nil,      nil}),tol);
                    
% unoEH=round(tt_blockwise({uno,     nil3;
%                         nil3,    uno}),tol);

% A     = round(sc_curlH*inv_mu_All/mu0*sc_curlE-omega^2*eps_All*eps0,tol);
% A_EH=round(tt_blockwise({1j*k0*eps_All,-sc_curlH; sc_curlE, 1j*k0*uno}),tol);

    if doUPML
        switch(units)
        case 'SI'
            A=curlH*inv_mu_All*curlE/mu0-omega^2*eps_All*eps0;
            A_EH=[1j*omega*eps0*eps_All,-curlH;...
            curlE, 1j*omega*mu0*mu_All];
        case 'Gauss'
            Aff=curlH*inv_mu_All*curlE-k0^2*eps_All;
            A_EH=[1j*k0*eps_All,-curlH;...
            sc_curlE, 1j*k0*mu0*uno];
        end
    else
        switch(units)
        case 'Gauss'
            A=sc_curlH*inv_mu_All*sc_curlE-k0^2*eps_All;
            A_EH=tt_blockwise({1j*k0*eps_All,-sc_curlH;...
            sc_curlE, 1j*k0*mu0*uno});
        case 'SI'
            A=sc_curlH*inv_mu_All*sc_curlE/mu0-omega^2*eps_All*eps0;
            A_EH=tt_blockwise({1j*omega*eps0*eps_All,-sc_curlH;...
            sc_curlE, 1j*omega*mu0*uno});
        end
    end


end

%% qtt-solver
if QTTsolverON
    
% fh_planewave=@(xyz) exp(-1j*(xyz(1)*k0x+xyz(2)*k0y));
% mask_amen=amen_cross(2*ones(d, 1),@(bin) fh_planewave(binarr2xyz(bin-1,...
%     dx,dy,dz,hx,hy,hz)), 1e-14,'nswp',40);

%generating rhs via ultrafast pre-defined tt_templates (much faster than
%amen_cross)
fexpX=round(tt_sin_cos(dx,k0x,1)+1j*tt_sin_cos(dx,k0x,0),tol);
fexpY=round(tt_sin_cos(dy,k0y,1)+1j*tt_sin_cos(dy,k0y,0),tol);
SOURCEFIELD1=mtkron(tt_ones(2,dz),fexpY,fexpX);

ttSOURCEFIELDE3_L=round(tt_blockwise({tt_zeros(2*ones(1,d));tt_zeros(2*ones(1,d))...
    ;SOURCEFIELD1}),tol);
ttSOURCEFIELDH3_L=round(tt_blockwise({tt_zeros(2*ones(1,d))...
    ;-SOURCEFIELD1;tt_zeros(2*ones(1,d))}),tol);

ttSOURCEFIELD3_L_EH=round(tt_blockwise({ttSOURCEFIELDE3_L;ttSOURCEFIELDH3_L}),tol);
b=round(round(round(ttQ*A,tol)-round(A*ttQ,tol),tol)*ttSOURCEFIELDE3_L,tol);

b_EH=round(round(round(ttQ_EH*A_EH,tol)-round(A_EH*ttQ_EH,tol),tol)*ttSOURCEFIELD3_L_EH,tol);

% bf3_L=full(ttb3_L);
% Pttb3_L = round(Ptt*tfsf_b3_L,tol);
% PA      = round(Ptt*A,tol);
% APPA    = round(PA'*PA,tol);
% APPb    = round(PA'*Pttb3_L,tol);
% 
% ttF3_L=round(amen_solve2(APPA,APPb,tol,'local_prec','rjacobi',...
%     'nswp',40,'resid_damp',2,'local_restart',60 ),tol);
% ttF3_L_2=round(amen_solve2(PA,Pttb3_L,1e-12,'local_prec','rjacobi',...
%     'nswp',100,'resid_damp',1,'local_restart',60 ),tol);
% ttF3_L_3=round(amen_solve2(A,tfsf_b3_L,1e-10,'local_prec','rjacobi',...
%     'nswp',100,'resid_damp',2,'local_restart',60 ),tol);

% % % AA_EH=round(A_EH'*A_EH,tol)
% % % Ab_EH=round(A_EH'*b_EH,tol)
% % % % sol.x={}
% % % % tol=1e-8;
% % % % %max_full_size 30000 is still smooth and performes well
% % % % 
% % % sol.x=round(amen_solve2(AA_EH,Ab_EH,tol,...
% % %         'local_prec','rjacobi','max_full_size',50000,'nswp',30,'resid_damp',1.3,...
% % %          'local_restart',5 ,'rmax',20),tol)
% % % sol.res=norm(AA_EH*sol.x-Ab_EH)/norm(Ab_EH)
% test_tensorcompression
QTTKrylov=0
if QTTKrylov
[U_gmres,td_gmres] = tt_gmres(core(A), core(b), tol, 1, 100, tol, tol, [], [], [], [], 1);
x_krylov=oldcell2core(U_gmres);
end

wavesol.x=amen_solve2(A,b,tol,'local_prec',...
     'rjacobi','max_full_size',100000,'nswp',15,...
     'resid_damp',2,'local_restart',10 ,'rmax',20);
 
wavesol.res=norm(A*wavesol.x-b)/norm(b)

% AA=round(A'*A,tol)
% Ab=round(A'*b,tol)
% wavesol_symm.x=amen_solve2(AA,Ab,tol,'local_prec',...
%      'rjacobi','max_full_size',100000,'nswp',30,...
%      'resid_damp',2,'local_restart',10 ,'rmax',20)
% wavesol_symm.res=norm(A*wavesol_symm.x-b)/norm(b)

% 
% field3L=full(wavesol.x);
% % 
% % % 
% % %  field3L_EH=full(sol.x);
%  field3D_EH=
%  (field3L_EH(2*N+1:3*N),params);
 
% % % 
% % %  QTT_EH_refl=abs(sol.x([xyz2bin([ xtest, ytest ,ztest],params)+1,...
% % %     CONST.AXIS.Z,CONST.FIELD.E]))
 
% norm(sol.x([xyz2bin([ xtest, ytest ,ztest],params)+1,...
%     CONST.AXIS.Z,CONST.FIELD.E])-field3D_EH(xtest+1,ytest+1,ztest+1))


%  F3L_visualize(field3L_EH,params,'field','Ez','fig',3,'name','QTT solver - EH');
 field3L=full(wavesol.x);
 field3D=F1_LtoF1_3D(field3L(2*N+1:3*N),params);
 
 QTT_wave_refl= abs(wavesol.x([xyz2bin([ xtest, ytest ,ztest],params)+1,...
    CONST.AXIS.Z])) 
%abs(field3D(floor(2^dx*tfsffrac)+2,2^(dy-1),2^(dz-1)))
 F3L_visualize(field3L,params,'field','Ez','fig',21,'name','QTT solver - wave equation');
end
%% plotting3d
% 
% letmePlot3d=0;
% if letmePlot3d
% figure()
% iso_mask=isosurface(X,Y,Z,fullmask1_3D);
% p=patch(iso_mask);
% p.FaceAlpha=0.9;
% p.EdgeColor   =   'none';
% daspect([1,1,1])
% view(3); axis tight
% p.FaceColor   =   'yellow';
% camlight 
% lighting gouraud
% end
%% converting geometry into FULL format and visualising via voxel method-----------------------------
% 
%  V=full(eps_Ex);
%  V=reshape(V,2^dz, 2^dy ,2^dx);
%  V=permute(V,[3 2 1]);
%  reshape
%  HA=vol3d('cdata',V,'texture','3D');
%  view(3);  
%  axis tight;  daspect([1 1 1])

