%By Alexey Boyko, Skolkovo Institute of Science and Technology,

%a.boyko@skoltech.ru
%alexey.boyko@skolkovotech.ru
%Moscow, Russia, 2015
%FDFD

%weirdly enough, only if one choses the tf\sf boundary at PML boundary,
%everything works

time.start=cputime;

clear all;
format long;
CONST.AXIS.X=1;
CONST.AXIS.Y=2;
CONST.AXIS.Z=3;
CONST.FIELD.E=1;
CONST.FIELD.H=2;
%% Setting up general constants -------------------------------
tol=1e-6;
SparseSolverON=1;
QTTsolverON=1;
%% Setting up the grid parameters------------------------------

dx=7;dy=7;dz=1;
d=dx+dy+dz;
hx=1;hy=1;hz=1;
h=min([hx hy hz]);

% Lx=1;Ly=1;Lz=1;
% hx=Lx/(2^dx); hy=Ly/(2^dy);hz=Lz/(2^dz);

doUPML=0 %UPML doesn't work for som

tfsffrac=0.75

XisPeriodic=0;
YisPeriodic=1;
ZisPeriodic=1;

%% Setting up physical envinronment and geometry and calculating in on qtt lattice-------
time.geom=cputime;
eps=4;
mu=1.0;

lambdas=(2^(dx-4)+2^(dx-8));
depth=ceil((2^(dx-3)+2^(dx-7))/2);%
% lambdas=depth;%.*[0.5:1.5/50:2];


params.tol=tol;
params.dx=dx;
params.dy=dy;
params.dz=dz;
params.d=params.dx+params.dy+params.dz;
params.hx=hx;
params.hy=hy;
params.hz=hz;
% params.Lx=2^params.dx-1;
% params.Ly=2^params.dy-1;
% params.Lz=2^params.dz-1;

% params_2X.tol=tol;
% params_2X.dx=dx+1;
% params_2X.dy=dy+1;
% params_2X.dz=dz+1;
% params_2X.d=params_2X.dx+params_2X.dy+params_2X.dz;
% params_2X.Lx=2^params_2X.dx-1;
% params_2X.Ly=2^params_2X.dy-1;
% params_2X.Lz=2^params_2X.dz-1;
% params_2X.hx=hx;
% params_2X.hy=hy;
% params_2X.hz=hz;

sf=@(x)sparse(full(x));

%eps=4
%8 7
% lambda=34
%holes=32

%eps=4
%9 8
% lambda=67.0
%holes = 64
%% Geometry
%both of index boundaries are included(!) in the box object 

% xb1_1=2^(dx-1)-floor(depth/2);xb1_2=2^(dx-1)+floor(depth/2);
% yb1_1=1;yb1_2=2^(dy);
% zb1_1=1;zb1_2=2^dz;
% box1=geom_box(dx,dy,dz,tol,xb1_1,xb1_2,yb1_1,yb1_2,zb1_1,zb1_2);

% hole=2*lambda;

holes=2^(dy-4);
geom={};

xb1=2^(dx-1)-floor(depth/2);
xb2=2^(dx-1)+floor(depth/2);

zb1=1;
zb2=2^dz;

yb1=1;%ceil(.1*2^dy);
yb2=2^dy;%ceil(.9*2^dy);

mask_of_coords=tt_zeros([2*ones(1,d)]);

nboxmax=(2^dy)/holes;

% for i=1:2:nboxmax
% geom{i}=geom_box(params,xb1,xb2,i*holes-floor(holes/2)+1,i*holes+floor(holes/2),zb1,zb2);
% mask_of_coords=round(mask_of_coords+geom{i},tol);
% end

geom{nboxmax+1}=purify(geom_box(params,xb1-depth,xb1-1,yb1,yb2,zb1,zb2));

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
xb3_1=1;xb3_2=2^dx*tfsffrac;
yb3_1=1;yb3_2=2^dy;
zb3_1=1;zb3_2=2^dz;


box_Q=geom_box(params,xb3_1,xb3_2,yb3_1,yb3_2,zb3_1,zb3_2); 
% % box_Q=geom_box(params_2X,xb3_1,xb3_2,yb3_1,yb3_2,zb3_1,zb3_2); %full 2X grid box
% % box_Q_Ex=tt_subtensor(box_Q,[1,dz+2,dz+dy+3],yee('Ex'));
% % box_Q_Ey=tt_subtensor(box_Q,[1,dz+2,dz+dy+3],yee('Ey'));
% % box_Q_Ez=tt_subtensor(box_Q,[1,dz+2,dz+dy+3],yee('Ez'));
% % box_Q_Hx=tt_subtensor(box_Q,[1,dz+2,dz+dy+3],yee('Hx'));
% % box_Q_Hy=tt_subtensor(box_Q,[1,dz+2,dz+dy+3],yee('Hy'));
% % box_Q_Hz=tt_subtensor(box_Q,[1,dz+2,dz+dy+3],yee('Hz'));

% box_Q_m=diag(box_Q);
% nil=ttm_zeros(2,dx+dy+dz);
% box_Q_All=round(tt_blockwise({box_Q_m,      nil,       nil ;
%                             nil,           box_Q_m,  nil;
%                             nil,           nil,       box_Q_m     }),tol)

if SparseSolverON 
    spbox_Q_Ex=sf(box_Q);
    spbox_Q_Ey=sf(box_Q);
    spbox_Q_Ez=sf(box_Q);

%     spbox_Q_Ex=sf(box_Q_Ex);
%     spbox_Q_Ey=sf(box_Q_Ey);
%     spbox_Q_Ez=sf(box_Q_Ez);
    Q=diag([spbox_Q_Ex;spbox_Q_Ey;spbox_Q_Ez]);
    
    spbox_Q_Hx=sf(box_Q);
    spbox_Q_Hy=sf(box_Q);
    spbox_Q_Hz=sf(box_Q);    
%     spbox_Q_Hx=sf(box_Q_Hx);
%     spbox_Q_Hy=sf(box_Q_Hy);
%     spbox_Q_Hz=sf(box_Q_Hz);
    Q_EH=diag([spbox_Q_Ex;spbox_Q_Ey;spbox_Q_Ez;spbox_Q_Hx;spbox_Q_Hy;spbox_Q_Hz]);
end
if QTTsolverON
    
    box_3Q_E=round(tt_blockwise({box_Q;box_Q;box_Q}),tol);
    box_3Q_H=round(tt_blockwise({box_Q;box_Q;box_Q}),tol);    
%     box_3Q_E=round(tt_blockwise({box_Q_Ex;box_Q_Ey;box_Q_Ez}),tol);
%     box_3Q_H=round(tt_blockwise({box_Q_Hx;box_Q_Hy;box_Q_Hz}),tol);
    box_6Q_EH=round(tt_blockwise({box_3Q_E;box_3Q_H}),tol);

    ttQ=diag(box_3Q_E);
    ttQ_EH=diag(box_6Q_EH);
end

%tt_ones(2,d)% round(box1+box2,tol);
%% running experiments
% this is a main loop for a series of experiments
for exp_i=1:numel(lambdas)
    lambda=lambdas(exp_i);
    tolx=(2^dx/lambda)^2;
    k0=2*pi/lambda;
    inc_angle=30/360*2*pi;
    inc_angle2=0/360*2*pi;

    k0x=k0*cos(inc_angle)*cos(inc_angle2)
    k0y=k0*sin(inc_angle)*cos(inc_angle2)
    k0z=k0*sin(inc_angle2)

    BlochBCS_EX=exp(-1j*k0x*(2^dx));
    BlochBCS_HX=BlochBCS_EX';
    BlochBCS_EY=exp(-1j*k0y*(2^dy));
    BlochBCS_HY=BlochBCS_EY';
    BlochBCS_EZ=exp(-1j*k0z*(2^dz));
    BlochBCS_HZ=BlochBCS_EY';
    
    

    params.k0=k0;

    pml_thickness=1*floor(lambda);

    pml.xlowwidth=pml_thickness;%%3*l ambda;
    pml.xhighwidth=pml_thickness;%3*lambda;
    pml.ylowwidth=pml_thickness;%3*lambda;
    pml.yhighwidth=pml_thickness;%3*lambda;
    pml.zlowwidth=0;
    pml.zhighwidth=0;
    pml.R=exp(-25);
    pml.m=2.5;

    %point for probing magnitude of reflected wave
     xtest=2^dx-pml.xhighwidth-2;
     xtest_trans=pml.xlowwidth+2;
     ytest=2^(dy-1);
     if dz>0
        ztest=2^(dz-1);
     else
        ztest=1
     end

    %% creating PML

    if doUPML == 1
        [syz_x,sxz_y,sxy_z,invsxz_y,invsyz_x,invsxy_z]=createUPML2(params,pml);
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


    %% sparse source creation
    time.sp.src=cputime;
    N=2^d;
    % planewave1_3D  =  @(x)exp(-1j*k0*x ).* double(x<200);
    % planewave1_3D  =  @(x)exp(-1j*k0*x ).* double(x>2*pml_xlow&(x<(2^dx-2*pml_xhigh-1)));

    if (d<=23)&&SparseSolverON
        planewave1_3D  =  @(x,y) exp(-1j*(k0x*x+k0y*y));
        %SOURCEFIELD1_3D=10*src_pointsource(k0,X,Y,Z,dx,dy,dz);
        SOURCEFIELD1_3D=planewave1_3D(X,Y);
        SOURCEFIELD3_L=[zeros(N,1);zeros(N,1);F1_3DtoF1_L(SOURCEFIELD1_3D)];
        SOURCEFIELD3_L_EH=[SOURCEFIELD3_L;zeros(N,1);-F1_3DtoF1_L(SOURCEFIELD1_3D);zeros(N,1)];
    end

    %% Assemblying Maxwell operator full\sparse format------------------------------
    time.sp.oper=cputime;
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
        if dz>0
            DEZ=sparse(diag(ones(2^dz-1,1),1)-diag(ones(2^dz,1),0));
            DEZ(2^dz,1)=ZisPeriodic*BlochBCS_EZ; %periodic 
        else
            DEZ=sparse(1);
        end
        %so right now X is the largest harmonic, Z is the highest-freq harmonic
        DHX=sparse(-diag(ones(2^dx-1,1),-1)+diag(ones(2^dx,1),0));
        DHX(1,2^dx) = -XisPeriodic*BlochBCS_HX;
        DHY=sparse(-diag(ones(2^dy-1,1),-1)+diag(ones(2^dy,1),0));
        DHY(1,2^dy) = -YisPeriodic*BlochBCS_HY;
        if dz>0
            DHZ=sparse(-diag(ones(2^dz-1,1),-1)+diag(ones(2^dz,1),0));
            DHZ(1,2^dz) = -ZisPeriodic*BlochBCS_HZ;
        else
            DHZ=sparse(1);
        end

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
            A_EH_ff=[1j*k0*eps_All_f,-curlH_f;...
            curlE_f, 1j*k0*speye(size(eps_All_f))];
        else
            Aff=sc_curlH_f*inv_mu_All_f*sc_curlE_f-k0^2*eps_All_f;
            A_EH_ff=[1j*k0*eps_All_f,-sc_curlH_f;...
            sc_curlE_f, 1j*k0*speye(size(eps_All_f))];
        end
        % A_empty_ff=curlH_f*curlE_f-k0^2*speye(2^d * 3) ;
        % A_empty_ff_scpml=sc_curlH_f*sc_curlE_f-k0^2*speye(2^d * 3) ;


        % 
        % A_empty_EH_ff=[1j*k0*speye(size(eps_All_f)), -sc_curlH_f;...
        %                 sc_curlE_f, 1j*k0*speye(size(eps_All_f))];



        % A_empty_ff_upml=curlH_f*coefs_inv_empty_pml*curlE_f-k0^2*coefs_empty_pml;

        % P=inv_eps_All_f;

        % Krylov solver for sparse
        % [L,U]   =   ilu(A)
        % src=zeros(2^dx,2^dy,2^dz);
        % src(2*pml_xlow,:,:)=1;
        %  bf3_L=[zeros(N,1);zeros(N,1);F1_3DtoF1_L(src)];
        %  bf3_L=Aff*SOURCEFIELD3_L;
        %  bf3_L_EH=A_EH_ff*SOURCEFIELD3_L_EH;

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
    time.sp.solve=cputime;
    F3_L=Aff\b_ff;

    F3_L_EH=A_EH_ff\b_EH_ff;
    time.sp.solved=cputime;

    ftt=tt_tensor( F3_L,tol,[2*ones(1,d),3])
     F3_L=full(ftt);
     F3_L_3d =    F1_LtoF1_3D(F3_L((2*N+1):(3*N)),params);

    fcross_refl = F3_L_3d(xtest,:,:);
    fcross_trans = F3_L_3d(xtest_trans,:,:);

    crossection_refl =reshape(fcross_refl,[numel(fcross_refl),1]);
    crossection_trans=reshape(fcross_trans,[numel(fcross_refl),1]);
    Rf(exp_i)=sum(abs(crossection_refl ).^2)/numel(crossection_refl)
    Tf(exp_i)=sum(abs(crossection_trans).^2)/numel(crossection_trans)
    Ef(exp_i)=Tf(exp_i)+Rf(exp_i)

    Fz_s_3d_EH = F1_LtoF1_3D(F3_L_EH((2*N+1):(3*N)),params);

    F3L_visualize(F3_L_EH,params,'field','Ez','fig',1,'name','Direct solver - EH');
    EH_refl=abs(Fz_s_3d_EH(xtest+1,ytest+1,ztest+1))

    F3L_visualize(F3_L,params,'field','Ez','name','Direct solver - Vector Wave','func','re');
    hold on
    contour(fullmask1_3D(:,:,2^(dz)),2,'color','White');
    hold off
    
    

    vectorwave_refl=abs(F3_L_3d(xtest,ytest,ztest))

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

    spsolution=time.sp.solved-time.sp.solve;

end
%% Assemblying Maxwell operator it qtt format------------------------------
if QTTsolverON
    time.qtt.oper=cputime;
    nil=ttm_zeros(2,dx+dy+dz);
    
    tt_eps_3_L=round((mask_of_coords*eps + (1-mask_of_coords)),tol)
    eps_Ex=tt_eps_3_L;
    eps_Ey=tt_eps_3_L;
    eps_Ez=tt_eps_3_L;
    
    if doUPML
        eps_Ex=round(eps_Ex.*syz_x,tol);
        eps_Ey=round(eps_Ey.*sxz_y,tol);
        eps_Ez=round(eps_Ez.*sxy_z,tol);
    end

    eps_Ex_m=diag(eps_Ex);
    eps_Ey_m=diag(eps_Ey);
    eps_Ez_m=diag(eps_Ez);

    eps_All=round(tt_blockwise({eps_Ex_m,      nil,       nil ;
                                nil,           eps_Ey_m,  nil;
                                nil,           nil,       eps_Ez_m     }),tol)
% keyboard
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

    if doUPML
        inv_mu_Ex=round(inv_mu_Ex.*invsyz_x,tol);
        inv_mu_Ey=round(inv_mu_Ey.*invsxz_y,tol);
        inv_mu_Ez=round(inv_mu_Ez.*invsxy_z,tol);
    end
    inv_mu_Ex_m  =  diag(inv_mu_Ex);
    inv_mu_Ey_m  =  diag(inv_mu_Ey);
    inv_mu_Ez_m  =  diag(inv_mu_Ez);

    inv_mu_All=round(tt_blockwise({inv_mu_Ex_m, nil,          nil;
                                   nil,         inv_mu_Ey_m,  nil;
                                   nil,         nil,          inv_mu_Ez_m     }),tol);             

    %warning! mtkron takes tensor product is reverse order
    
    ttDEX=1/hx * fd_tt_DE(2,dx,XisPeriodic,tol,BlochBCS_EX);
    ttDEY=1/hy * fd_tt_DE(2,dy,YisPeriodic,tol,BlochBCS_EY);
    if dz>0
        ttDEZ=1/hz * fd_tt_DE(2,dz,ZisPeriodic,tol,BlochBCS_EZ);
    else
        ttDEZ=tt_matrix(1);
    end
    deltaEx=purify(mtkron(tt_eye(2,dz),tt_eye(2,dy),ttDEX));
    deltaEy=purify(mtkron(tt_eye(2,dz),ttDEY,tt_eye(2,dx)));
    deltaEz=purify(mtkron(ttDEZ,tt_eye(2,dy),tt_eye(2,dx)));

    %so right now X is the largest harmonic, Z is the highest-freq harmonic
    
    ttDHX=1/hx * fd_tt_DH(2,dx,XisPeriodic,tol,BlochBCS_HX);
    ttDHY=1/hy * fd_tt_DH(2,dy,YisPeriodic,tol,BlochBCS_HY);
    if dz>0
        ttDHZ=1/hz * fd_tt_DH(2,dz,ZisPeriodic,tol,BlochBCS_HZ);
    else
        ttDHZ=tt_matrix(1);
    end       
    deltaHx=purify(mtkron(tt_eye(2,dz),tt_eye(2,dy),ttDHX));
    deltaHy=purify(mtkron(tt_eye(2,dz),ttDHY,tt_eye(2,dx)));
    deltaHz=purify(mtkron(ttDHZ,tt_eye(2,dy),tt_eye(2,dx)));

    uno=round(tt_blockwise({tt_eye(2,d),     nil,      nil;
                            nil,     tt_eye(2,d),      nil;
                            nil,     nil,      tt_eye(2,d)}),tol);
    nil3=round(tt_blockwise({nil,     nil,      nil;
                            nil,     nil,      nil;
                            nil,     nil,      nil}),tol);

    if ~doUPML
        sc_deltaEx=diag(invsx)*deltaEx;
        sc_deltaEy=diag(invsy)*deltaEy;
        sc_deltaEz=diag(invsz)*deltaEz;

        sc_deltaHx=diag(invsx)*deltaHx;
        sc_deltaHy=diag(invsy)*deltaHy;
        sc_deltaHz=diag(invsz)*deltaHz;
        sc_curlE=round(tt_blockwise({nil,           -sc_deltaEz,       sc_deltaEy ;
                        sc_deltaEz,        nil,            -sc_deltaEx;
                        -sc_deltaEy,       sc_deltaEx,        nil     }),tol);

        sc_curlH=round(tt_blockwise({nil,           -sc_deltaHz,       sc_deltaHy ;
                            sc_deltaHz,        nil,            -sc_deltaHx;
                            -sc_deltaHy,       sc_deltaHx,        nil     }),tol);  

        A     = round(sc_curlH*inv_mu_All*sc_curlE-k0^2*eps_All,tol);


        % unoEH=round(tt_blockwise({uno,     nil3;
        %                         nil3,    uno}),tol);


        A_EH=round(tt_blockwise({1j*k0*eps_All,-sc_curlH; sc_curlE, 1j*k0*uno}),tol);
    else
        % compressing deltas doesnt help much


        curlE=tt_blockwise({nil,           -deltaEz,       deltaEy ;
                            deltaEz,        nil,            -deltaEx;
                            -deltaEy,       deltaEx,        nil     });

        curlH=tt_blockwise({nil,           -deltaHz,       deltaHy ;
                            deltaHz,        nil,            -deltaHx;
                            -deltaHy,       deltaHx,        nil     });   
        
        A     = round(curlH*inv_mu_All*curlE-k0^2*eps_All,tol);
        A_EH=round(tt_blockwise({1j*k0*eps_All,-curlH; curlE, 1j*k0*uno}),tol);
    end              

    % A= round(curlH*curlE-k0^2*diag(tt_ones([2*ones(1,d) 3])),tol);

    % A_empty = round(curlH*curlE-k0^2*diag(tt_ones([2*ones(1,d) 3])),tol);

    % A_empty_scpml = round(sc_curlH*sc_curlE-k0^2*diag(tt_ones([2*ones(1,d) 3])),tol);


    % A_EH_cor=round(A_EH+unoEH,tol)

end

%% qtt-solver sources
    if QTTsolverON
        time.qtt.src=cputime;    
        % fh_planewave=@(xyz) exp(-1j*(xyz(1)*k0x+xyz(2)*k0y));
        % mask_amen=amen_cross(2*ones(d, 1),@(bin) fh_planewave(binarr2xyz(bin-1,...
        %     dx,dy,dz,hx,hy,hz)), 1e-14,'nswp',40);

        %generating rhs via ultrafast pre-defined tt_templates (much faster than
        %amen_cross)
        fexpX=round(tt_sin_cos(dx,k0x,1)-1j*tt_sin_cos(dx,k0x,0),tol);
        fexpY=round(tt_sin_cos(dy,k0y,1)-1j*tt_sin_cos(dy,k0y,0),tol);
        if dz>1
            fexpZ=round(tt_sin_cos(dz,k0z,1)-1j*tt_sin_cos(dz,k0z,0),tol);
        else
            fexpZ=tt_ones(2,dz)
        end
        SOURCEFIELD1=purify(mtkron(fexpZ,fexpY,fexpX));

        ttSOURCEFIELDE3_L=round(tt_blockwise({tt_zeros(2*ones(1,d));tt_zeros(2*ones(1,d))...
            ;SOURCEFIELD1}),tol);
        ttSOURCEFIELDH3_L=round(tt_blockwise({tt_zeros(2*ones(1,d))...
            ;-SOURCEFIELD1;tt_zeros(2*ones(1,d))}),tol);
        % ttSOURCEFIELD3_L=tt_tensor(SOURCEFIELD3_L,tol,[2*ones(1,d),3]);
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
        

%     AA_EH=round(A_EH'*A_EH,tol)
%     Ab_EH=round(A_EH'*b_EH,tol)
%     sol.x=amen_solve2(AA_EH,Ab_EH,tol,'local_prec',...
%         'rjacobi','max_full_size',MFS,'nswp',30,...
%         'resid_damp',1.5,'local_restart',1 ,'rmax',70,'trunc_norm','fro','x0',Ab_EH);
% 
%        
        % % % % %max_full_size 30000 is still smooth and performes well
        
        %% qtt solution
  
%         sol.x=round(amen_solve2(AA_EH,Ab_EH,tol,...
%                 'local_prec','rjacobi','max_full_size',10000,'nswp',20,'resid_damp',1.3,...
%                  'local_restart',5 ,'rmax',30),tol)
%         sol.res=norm(AA_EH*sol.x-Ab_EH)/norm(Ab_EH)
%         sol.x=round(amen_solve2(AA_EH,Ab_EH,tol,...
%                 'local_prec','rjacobi','max_full_size',10000,'nswp',7,'resid_damp',1.3,...
%                  'local_restart',5 ,'rmax',40,'x0',),tol)
%         sol.res=norm(AA_EH*sol.x-Ab_EH)/norm(Ab_EH)        

        
        % test_tensorcompression
        QTTKrylov=0
        if QTTKrylov
            [U_gmres,td_gmres] = tt_gmres(core(A), core(b), tol, 1, 100, tol, tol, [], [], [], [], 1);
            x_krylov=oldcell2core(U_gmres);
        end

        amentols=10.^-[1:14]
        time.qtt.solve=cputime; 
        
        tol
        
        MFS=15000;
        startrank=20;
        wavesol_exp.x{1}=amen_solve2(A,b,tol,'local_prec',...
        'rjacobi','max_full_size',MFS,'nswp',14,...
        'resid_damp',1.5,'local_restart',1 ,'rmax',startrank,'trunc_norm','fro','x0',b);
        time.qtt.solved{1}=cputime;
   %     wavesol_exp.res{1}=norm(A*wavesol_exp.x{1}-b)/norm(b)
        for i=2:5
            wavesol_exp.x{i}=amen_solve2(A,b,tol,'local_prec',...
             'rjacobi','max_full_size',MFS,'nswp',6,...
             'resid_damp',1.7,'kickrank',6,'local_restart',1 ,'rmax',startrank+10*(i-1),'trunc_norm','fro','x0',wavesol_exp.x{i-1});
            time.qtt.solved{i}=cputime;
      %      wavesol_exp.res{i}=norm(A*wavesol_exp.x{i}-b)/norm(b)
        end
            
        
        
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


        % % % 
        % % %  QTT_EH_refl=abs(sol.x([xyz2bin([ xtest, ytest ,ztest],params)+1,...
        % % %     CONST.AXIS.Z,CONST.FIELD.E]))

        % norm(sol.x([xyz2bin([ xtest, ytest ,ztest],params)+1,...
        %     CONST.AXIS.Z,CONST.FIELD.E])-field3D_EH(xtest+1,ytest+1,ztest+1))



        %  F3L_visualize(field3L_EH,params,'field','Ez','fig',3,'name','QTT solver - EH');

        % % %  
        % % %  F3L_visualize(field3L_EH,params,'field','Ez','fig',3);

        wavesol.x=wavesol_exp.x{i}
        %  field3L=full(wavesol.x);
        %  field3D=F1_LtoF1_3D(field3L(2*N+1:3*N),params);
        QTT_wave_refl= abs(wavesol.x([xyz2bin([ xtest, ytest ,ztest],params)+1,...
            CONST.AXIS.Z])) 
        % %abs(field3D(floor(2^dx*tfsffrac)+2,2^(dy-1),2^(dz-1)))
        %  F3L_visualize(field3L,params,'field','Ez','fig',4,'name','QTT solver - wave equation','func','im');
        % test_subfield2
        wz=tt_subtensor(wavesol.x,d+1,CONST.AXIS.Z);

        ttv_crossection_refl=tt_subtensor(wz,dz+dy+1:dz+dy+dx ,[1+de2bi(xtest,dx)]);
        crossection_refl=full(ttv_crossection_refl);

        ttv_crossection_trans=tt_subtensor(wz,dz+dy+1:dz+dy+dx ,[1+de2bi(xtest_trans,dx)]);
        crossection_trans=full(ttv_crossection_trans);

        fcross_refl =reshape(crossection_refl,[2^dz 2^dy]);
        fcross_trans=reshape(crossection_trans,[2^dz 2^dy]);

%         log_err=log10(abs(1-E))
% 
%         qttsolution=time.qtt.solved-time.qtt.solve;
        
        RR(exp_i)=sum(abs(crossection_refl ).^2)/numel(crossection_refl)
        TT(exp_i)=sum(abs(crossection_trans).^2)/numel(crossection_trans)
        EE(exp_i)=TT(exp_i)+RR(exp_i);
        test_subfield2
    end



end

% load('../')
% figure;
% plot(lambdas/depth,R,'x-r');
% hold on;
% plot(lambdas/depth,T,'x-b');
% plot(lambdas/depth,RR,'o--r');
% plot(lambdas/depth,TT,'o--b');
% hold off
% legend('R','T','Rqtt','Tqtt')

if exist('spsolution')
spsolution
end

if exist('qttsolution')
qttsolution
end
%% plotting3d

letmePlot3d=0;
if letmePlot3d
figure()
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
%% converting geometry into FULL format and visualising -----------------------------
% 
%  V=full(eps_Ex);
%  V=reshape(V,2^dz, 2^dy ,2^dx);
%  V=permute(V,[3 2 1]);
%  reshape
%  HA=vol3d('cdata',V,'texture','3D');
%  view(3);  
%  axis tight;  daspect([1 1 1])

