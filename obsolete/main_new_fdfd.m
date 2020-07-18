%By Alexey Boyko, Skolkovo Institute of Science and Technology,
%a.boyko@skoltech.ru
%alexey.boyko@skolkovotech.ru
%Moscow, Russia, 2015
%FDFD 
% 

%weirdly enough, only if one choses the tf\sf boundary at PML boundary,
%everything works
clear all;

%% Setting up general constants -------------------------------
tol=1e-6;

SparseSolverON=1;
QTTsolverON=0;
%% Setting up the grid parameters------------------------------
dx=8;dy=7;dz=1;
d=dx+dy+dz;

XisPeriodic=0;
YisPeriodic=1;
ZisPeriodic=1;


%% Setting up physical envinronment and geometry and calculating in on qtt lattice-------
eps=-16.0;
mu=1.0;
lambda=2^(dy-4);
k0=2*pi/lambda;

inc_angle=0/360*2*pi;

k0x=k0*cos(inc_angle);
k0y=k0*sin(inc_angle);

BlochBCS_EX=exp(-1j*k0x*(2^dx));
BlochBCS_HX=BlochBCS_EX';
BlochBCS_EY=exp(-1j*k0y*(2^dy));
BlochBCS_HY=BlochBCS_EY';

% Lx=1;Ly=1;Lz=1;
% hx=Lx/(2^dx); hy=Ly/(2^dy);hz=Lz/(2^dz);
hx=1;hy=1;hz=1;
h=min([hx hy hz]);

hole=2*lambda;
depth=2*lambda;

pml_thickness=20;
pml_xlow=pml_thickness;%%3*lambda;
pml_xhigh=pml_thickness;%3*lambda;
pml_ylow=pml_thickness;%3*lambda;
pml_yhigh=pml_thickness;%3*lambda;
pml_zlow=0;
pml_zhigh=0;
doUPML=0;
if doUPML == 1
[syz_x,sxz_y,sxy_z,invsxz_y,invsyz_x,invsxy_z]=createUPML(k0,...
    [pml_xlow,pml_xhigh,pml_ylow,pml_yhigh,pml_zlow,pml_zhigh],...
    tol,dx,dy,dz,hx,hy,hz);
else
      [invsx,invsy,invsz]=createSCPML(k0,...
        [pml_xlow,pml_xhigh,pml_ylow,pml_yhigh,pml_zlow,pml_zhigh],...
        tol,dx,dy,dz,hx,hy,hz);
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

xb1_1=2^(dx-1)-floor(depth/2);xb1_2=2^(dx-1)+floor(depth/2);
yb1_1=1;yb1_2=2^(dy-1)-2*floor(hole);
zb1_1=1;zb1_2=2^dz;
box1=geom_box(dx,dy,dz,tol,xb1_1,xb1_2,yb1_1,yb1_2,zb1_1,zb1_2);

xb2_1=2^(dx-1)-floor(depth/2);xb2_2=2^(dx-1)+floor(depth/2);
yb2_1=2^(dy-1)+-1*floor(hole);yb2_2=2^(dy);
zb2_1=1;zb2_2=2^dz;
box2=geom_box(dx,dy,dz,tol,xb2_1,xb2_2,yb2_1,yb2_2,zb2_1,zb2_2);
% 
% xb1_1=2^(dx-1)-floor(depth/2);xb1_2=2^(dx-1)+floor(depth/2);
% yb1_1=1;yb1_2=2^(dy);
% zb1_1=1;zb1_2=2^dz;
% box=geom_box(dx,dy,dz,tol,xb1_1,xb1_2,yb1_1,yb1_2,zb1_1,zb1_2);
% 
xb3_1=2^dx-60;xb3_2=2^dx;
yb3_1=1;yb3_2=2^dy;
zb3_1=1;zb3_2=2^dz;

 box_Q=geom_box(dx,dy,dz,tol,xb3_1,xb3_2,yb3_1,yb3_2,zb3_1,zb3_2);
% box_Q_m=diag(box_Q);
% nil=ttm_zeros(2,dx+dy+dz);
% box_Q_All=round(tt_blockwise({box_Q_m,      nil,       nil ;
%                             nil,           box_Q_m,  nil;
%                             nil,           nil,       box_Q_m     }),tol)
 spbox_Q=sparse(full(box_Q));
 
%  spbox_Q=speye(size(spbox_Q))-spbox_Q
 Q=diag([spbox_Q;spbox_Q;spbox_Q]);

mask_of_coords  =  box1+box2;%tt_ones(2,d)% round(box1+box2,tol);
%% source creation

[fullmask1_3D, X, Y, Z]  =  F1_LtoF1_3D(full(mask_of_coords), dx, dy, dz, hx, hy, hz);

N=2^d;
% planewave1_3D  =  @(x)exp(-1j*k0*x ).* double(x<200);
% planewave1_3D  =  @(x)exp(-1j*k0*x ).* double(x>2*pml_xlow&(x<(2^dx-2*pml_xhigh-1)));
 planewave1_3D  =  @(x,y) exp(-1j*(k0x*x+k0y*y));
% SOURCEFIELD1_3D=10*src_pointsource(k0,X,Y,Z,dx,dy,dz);
% 
SOURCEFIELD1_3D=planewave1_3D(X,Y);

SOURCEFIELD3_L=[zeros(N,1);zeros(N,1);F1_3DtoF1_L(SOURCEFIELD1_3D)];
%% Assemblying Maxwell operator full\sparse format------------------------------
if SparseSolverON
nil_f=sparse(2^d,2^d);

mask_of_coords_f=sparse(round(full(mask_of_coords)));

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
    Aff=curlH_f*inv_mu_All_f*curlE_f-k0^2*eps_All_f;
else
    Aff=sc_curlH_f*inv_mu_All_f*sc_curlE_f-k0^2*eps_All_f;
end
A_empty_ff=curlH_f*curlE_f-k0^2*speye(2^d * 3) ;
A_empty_ff_scpml=sc_curlH_f*sc_curlE_f-k0^2*speye(2^d * 3) ;
% A_empty_ff_upml=curlH_f*coefs_inv_empty_pml*curlE_f-k0^2*coefs_empty_pml;

% P=inv_eps_All_f;

% Krylov solver for sparse
% [L,U]   =   ilu(A)
% src=zeros(2^dx,2^dy,2^dz);
% src(2*pml_xlow,:,:)=1;
%  bf3_L=[zeros(N,1);zeros(N,1);F1_3DtoF1_L(src)];
bf3_L=A_empty_ff_scpml*1e4*SOURCEFIELD3_L;
%bf3_L(173780:196608)=0;

% keyboard
tfsf_bf3_L=(Q*Aff-Aff*Q)*bf3_L;
% [F3_L, flag, relres]=bicgstab(Aff, bf3_L, 1e-12,1000,[],[],SOURCEFIELD3_L); 
% flag
% bicgstab_res=norm(Aff*F3_L-bf3_L,1)

% Direct solver for sparse
% F3_L_source_pml=A_empty_ff_pml\bf3_L;
% bf3_L_1=A_empty_ff_pml
F3_L_orig=Aff\bf3_L;
F3_L=Aff\tfsf_bf3_L;
% dir_res=norm(Aff*F3_L1-bf3_L,1)
% 
% figKr=plot3fields(F3_L,bf3_L,2^(dz-1), dx, dy, dz, hx, hy, hz,1,'Krylov Solver');

N=2^(dx+dy+dz);

Fz_s_3d = F1_LtoF1_3D(F3_L((2*N+1):(3*N)),dx,dy,dz,hx,hy,hz);
figure();imagesc(real(Fz_s_3d(:,:,2^(dz-1))));
title('Ez')
axis equal tight
colorbar;

% Fz_s_3d_orig = F1_LtoF1_3D(F3_L_orig((2*N+1):(3*N)),dx,dy,dz,hx,hy,hz);
% figure();imagesc(real(Fz_s_3d_orig(:,:,2^(dz-1))));
% title('Ez')
% axis equal tight
% colorbar;

 figDir=plot3fields(F3_L,bf3_L,fullmask1_3D, 2^(dz-1), dx, dy, dz, hx, hy, hz,2,'Direct Solver')

end
%% Assemblying Maxwell operator it qtt format------------------------------
if QTTsolverON
nil=ttm_zeros(2,dx+dy+dz);


eps_Ex=round((mask_of_coords*eps + (1-mask_of_coords)).*syz_x,tol);
eps_Ey=round((mask_of_coords*eps + (1-mask_of_coords)).*sxz_y,tol);
eps_Ez=round((mask_of_coords*eps + (1-mask_of_coords)).*sxy_z,tol);

eps_Ex_m=diag(eps_Ex);
eps_Ey_m=diag(eps_Ey);
eps_Ez_m=diag(eps_Ez);

eps_All=round(tt_blockwise({eps_Ex_m,      nil,       nil ;
                            nil,           eps_Ey_m,  nil;
                            nil,           nil,       eps_Ez_m     }),tol)

%--for preconditioner----------------------------
                        
inv_eps_Ex=round(mask_of_coords*eps^(-1) + (1-mask_of_coords),tol);
inv_eps_Ey=round(mask_of_coords*eps^(-1) + (1-mask_of_coords),tol);
inv_eps_Ez=round(mask_of_coords*eps^(-1) + (1-mask_of_coords),tol);

inv_eps_Ex_m  =  diag(inv_eps_Ex);
inv_eps_Ey_m  =  diag(inv_eps_Ey);
inv_eps_Ez_m  =  diag(inv_eps_Ez);

inv_eps_All=round(tt_blockwise({inv_eps_Ex_m, nil,          nil;
                               nil,         inv_eps_Ey_m,  nil;
                               nil,         nil,          inv_eps_Ez_m     }),tol)           

%------------------------------------------------
                        
inv_mu_Ex=round((mask_of_coords*eps + (1-mask_of_coords)).*invsyz_x,tol);
inv_mu_Ey=round((mask_of_coords*eps + (1-mask_of_coords)).*invsxz_y,tol);
inv_mu_Ez=round((mask_of_coords*eps + (1-mask_of_coords)).*invsxy_z,tol);

inv_mu_Ex_m  =  diag(inv_mu_Ex);
inv_mu_Ey_m  =  diag(inv_mu_Ey);
inv_mu_Ez_m  =  diag(inv_mu_Ez);

inv_mu_All=round(tt_blockwise({inv_mu_Ex_m, nil,          nil;
                               nil,         inv_mu_Ey_m,  nil;
                               nil,         nil,          inv_mu_Ez_m     }),tol);             

%warning! mtkron takes tensor product is reverse order
deltaEx=mtkron(tt_eye(2,dz),                tt_eye(2,dy),                   1/hx * fd_tt_DE(2,dx,XisPeriodic,tol));
deltaEy=mtkron(tt_eye(2,dz),                1/hy * fd_tt_DE(2,dy,YisPeriodic,tol),    tt_eye(2,dx));
deltaEz=mtkron(1/hz * fd_tt_DE(2,dz,ZisPeriodic,tol), tt_eye(2,dy),                   tt_eye(2,dx));

%so right now X is the largest harmonic, Z is the highest-freq harmonic

deltaHx=mtkron(tt_eye(2,dz),                tt_eye(2,dy),               1/hx * fd_tt_DH(2,dx,XisPeriodic,tol));
deltaHy=mtkron(tt_eye(2,dz),                1/hy * fd_tt_DH(2,dy,YisPeriodic,tol),tt_eye(2,dx));
deltaHz=mtkron(1/hz * fd_tt_DH(2,dz,ZisPeriodic,tol), tt_eye(2,dy),               tt_eye(2,dx));

% compressing deltas doesnt help much

curlE=tt_blockwise({nil,           -deltaEz,       deltaEy ;
                    deltaEz,        nil,            -deltaEx;
                    -deltaEy,       deltaEx,        nil     });

curlH=tt_blockwise({nil,           -deltaHz,       deltaHy ;
                    deltaHz,        nil,            -deltaHx;
                    -deltaHy,       deltaHx,        nil     });   

% nil_nabla=round(tt_blockwise({nil,nil,nil;
%                               nil,nil,nil;
%                               nil,nil,nil}),tol); 

% Jx=tt_zeros(2,d);
% Jy=tt_zeros(2,d);
% Jz=source;

% b=-1j*round(diag(tt_blockwise({diag(Jx),nil,nil;
%                     nil,diag(Jy),nil;
%                     nil,nil,diag(Jz)})),tol)*k0; 

                          
A      =round(curlH*inv_mu_All*curlE-k0^2*eps_All,tol)
A_empty=round(curlH*curlE-k0^2*diag(tt_ones([2*ones(1,d) 3])),tol)
Ptt=inv_eps_All
% 
% Af=full(A);
% A_emptyf=full(A_empty);
% bf=full(b);
% Efield_f=Af\bf;

% Efield_f_r=reshape(Efield_f,[2^dz 2^dy 2^dx 3]);
% Efield_f_r=permute(Efield_f_r,[3 2 1 4]);

% 
% Efield=amen_solve2(A,b,tol,'nswp',100);
% [Efull,~,~,~]=QTTFIELDto3DFIELD( Efield(:,:,:,:,:,:,:,:,:,1),dx,dy,dz,hx,hy,hz );
end


%% qtt-solver
if QTTsolverON
ttSOURCEFIELD3_L=tt_tensor(SOURCEFIELD3_L,tol,[2*ones(1,d),3]);
ttb3_L=round(A_empty*ttSOURCEFIELD3_L,tol);
bf3_L=full(ttb3_L);
% Pttb3_L=round(Ptt*ttb3_L,tol);
% PA=round(Ptt*A,tol)
% % APPA=round(PA'*PA,tol)
% APPb=round(PA'*Pttb3_L,tol);
ttF3_L=amen_solve2(A,ttb3_L,tol,'local_prec','rjacobi','nswp',40 )
ttF3_L=round(ttF3_L,tol);
F3_L2=full(ttF3_L);
figQTT=plot3fields(F3_L2,bf3_L,2^(dz-1), dx, dy, dz, hx, hy, hz,3,'QTT Solver');
end
%% plotting

letmePlot3d=0;
if letmePlot3d

figure(1)
iso_mask=isosurface(X,Y,Z,fullmask1_3D);
p=patch(iso_mask);
p.FaceAlpha=0.3;
p.EdgeColor   =   'none';
daspect([1,1,1])
view(3); axis tight
p.FaceColor   =   'yellow';
camlight 
lighting gouraud
end
%% converting geometry into FULL format and visualising -----------------------------
% 
%  V=full(eps_Ex);
%  V=reshape(V,2^dz, 2^dy ,2^dx);
%  V=permute(V,[3 2 1]);
%  reshape
%  HA=vol3d('cdata',V,'texture','3D');
%  view(3);  
%  axis tight;  daspect([1 1 1])

