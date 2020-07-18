%this piece of code is intended to run qtt and fft and full-based Helmholtz
%2D equation, and compare solutions obtained by all methods

%case where Nx != Ny is not tested and probably will not work

%%Setting up environment constants
clear all;
%   close all;
experiment.dx0=9;
experiment.dxmax=experiment.dx0;

experiment.dxes=experiment.dx0:experiment.dxmax;
imax=experiment.dxmax-experiment.dx0+1;

%% Control block
DrawMasks=1;

FFTSolverOn=0;
DrawFFTSolutions = FFTSolverOn;
LOAD_FFT_ENTITY=0;

QTTSolverOn=1;
DrawQTTSolutions = 1;
CreateEntity=0;

WavelengthPerDomain=[15];

LOAD_MIE=0;
CALC_MIE=0;
DrawMIESolutions=0;
letmePlot3d=0;

if LOAD_FFT_ENTITY
	colloc_hardcode=matfile('datafiles\Gfft_5_12.mat');
end

if LOAD_MIE
    miefield_hardcode=matfile('datafiles\mie\data_mie_5_11.mat');
    data_mie=miefield_hardcode.data_mie;
end


%% Looping through experiment space
for i=1:numel(WavelengthPerDomain) 
%% 
params.grid.dx=experiment.dx0;
params.tol=1e-4;
%eee=@(ss) .1*(2.^(-8 - 2*(ss - 8)));
%tol=eee(dx)  %Vangovanie of tol
%% gridder block
params.grid.dy=params.grid.dx;
params.grid.dz=0;
params.grid.d=params.grid.dx+params.grid.dy+params.grid.dz;

params.grid.XLEFT=-.5; %use -2lambda, +2 lambda for comparison with Mie
params.grid.XRIGHT=.5;
params.grid.YLEFT=-.5;
params.grid.YRIGHT=.5;

params.grid.Lx=params.grid.XRIGHT-params.grid.XLEFT; 
params.grid.Ly=params.grid.YRIGHT-params.grid.YLEFT;
params.grid.Lz=1;

% |[ . ][ . ][ . ][ . ]|
% |<--------L--------->|

params.grid.Nx=2^params.grid.dx; 
params.grid.Ny=2^params.grid.dy;
params.grid.Nz=2^params.grid.dz;

params.grid.D=2+double(params.grid.dz>0); %dimensionality

params.grid.hx=params.grid.Lx/(params.grid.Nx);
params.grid.hy=params.grid.Ly/(params.grid.Ny);
params.grid.hz=params.grid.Lz/(params.grid.Nz);

%% RHS block
% params.rhs.lambda=params.grid.Lx / 40 * 361;
params.rhs.lambda=params.grid.Lx/WavelengthPerDomain(i);

% R=params.grid.Lx*0.25; %radius of a cylinder

% ks=10.^(0:.2:2)A()
params.rhs.k0=2*pi/params.rhs.lambda;
params.rhs.alpha=0/360*2*pi;

params.rhs.k0x=params.rhs.k0*cos(params.rhs.alpha);
params.rhs.k0y=params.rhs.k0*sin(params.rhs.alpha);

% % % chi=k0^2*(e_r-1);
%% technical block

params.quadratures.maxerrorjumps=20;
params.quadratures.softeningcoef=1/2;
if params.grid.d>22 && FFTSolverOn 
    error('trying to run FFT solver with d>22, aborting')
end

%% Create TT Entity (weirdly enough, needed mostly for FFT solver)
if ~LOAD_FFT_ENTITY && FFTSolverOn || QTTSolverOn ||CreateEntity

    if params.grid.D==2 
        ttGXY=abinitioG_entity(params);
        
%         GXY=reshape(full(ttGXY),[2*params.grid.Nx,2*params.grid.Ny]);
%         figure('name','Greens function over delta - QTT entity representation');imagesc(real(GXY));  colorbar;axis equal tight;    
     elseif params.grid.D==3
        Gh_xyz=@(deltax,deltay,deltaz) -1/6*(...
              exp(1j*k0*norm([deltax,deltay,deltaz] -[0 0 hz/2]))   / ( 4*pi*( norm([deltax,deltay,deltaz]-[0 0 hz/2]))) );
        Gh_r3=@(r) -1/6*(...
              exp(1j*k0*norm(r -[0 0 hz/2]))   / ( 4*pi*( norm(r-[0 0 hz/2]))) );
    %     ttGXYZ= amen_cross(2*ones(d+3, 1),@(bin) Gh_r3(bin2deltaxyz(bin-1,params)), tol,'nswp',40,'kickrank',12);  
    %     GXYZ=reshape(full(ttGXYZ),[2*Nx,2*Ny,2*Nz]);
          [Xdeltas,Ydeltas,Zdeltas]=meshgrid( (-2^dx:2^dx-1)*hx,...
            (-2^dy:2^dy-1)*hy,(-2^dz:2^dz-1)*hz);
         GXYZ=arrayfun(Gh_xyz,Xdeltas,Ydeltas,Zdeltas);
         ttGXYZ=tt_tensor(reshape(GXYZ,2*ones(1,d+3)),tol*1e-2)
    end
    
    
%     expdata{params.grid.dx}.WavelengthPerDomain=WavelengthPerDomain;
%     expdata{params.grid.dx}.ttGXYranks=ttGXY.r;
%     expdata{params.grid.dx}.tol=tol;
%     save(strcat(['ttGXY_expdata, wave=',num2str(WavelengthPerDomain),'.mat']),'','-v7.3')

end
%% Create or Load FFT Entity from TT Entity
if FFTSolverOn
    if params.grid.D==2
        if ~LOAD_FFT_ENTITY
            Gdeltas_fft=reshape(full(ttGXY),[2*params.grid.Nx,2*params.grid.Ny]);
            Gdeltas_fft=Gdeltas_fft(2:size(Gdeltas_fft,1),2:size(Gdeltas_fft,2)); %removing garbage padding
            Gdeltas_fft=circshift(circshift(Gdeltas_fft,-params.grid.Nx+1,1),-params.grid.Ny+1,2);
%            figure('name','Greens function over delta - QTT entity representation');imagesc(imag(GXY));  colorbar;axis equal tight;      
           figure('name','Greens function over delta - FFT entity representation');imagesc(real(Gdeltas_fft));  colorbar;axis equal tight;
        else
    % % % % % % % load 
            buf=colloc_hardcode.Gdeltas_fft_all(dx-5+1,1);
            Gdeltas_fft2=buf{1}*(hx*hy*2^d); % TEMP HX HY
%             clear buf;
            Gdeltas_fft=[rot90(Gdeltas_fft2,2),rot90(Gdeltas_fft2,1);rot90(Gdeltas_fft2,3),Gdeltas_fft2];        
            Gdeltas_fft(Nx,:)=[];
            Gdeltas_fft(:,Ny)=[];
            Gdeltas_fft=circshift(circshift(Gdeltas_fft,-Nx+1,1),-Ny+1,2);
%             clear Gdeltas_fft2;
    % % % % % % %         figure('name','Collocation integral from Mathematica - FFT entity');imagesc(real(Gdeltas_fft));  colorbar;axis equal tight;
        end
    elseif params.grid.D==3
        Gdeltas_fft=GXYZ(2:size(GXYZ,1),2:size(GXYZ,2),2:size(GXYZ,3));%removing garbage padding
        Gdeltas_fft=circshift(circshift(circshift(Gdeltas_fft,-Nx+1,1),-Ny+1,2),-Nz+1,3); 
    end
end
% % % % G_fft=blocktoeplitz(Gdeltas_fft);
% % % % norm(ffG-G_fft)
%% params.parts.eps constant
params.parts.eps=2.25 ;%-1-.29j%Ag_experimental(.361);
%-1.07-0.29j; % negative imaginary part for decay
%% Material mask forming 

% disp('Forming material from full format')

%___Mesh___


% [XX,YY]=meshgrid(XLEFT+ hx/2+ (0:(2^dx-1))*hx,YLEFT+ hy/2+(0:(2^dy-1))*hy);
INeedCircle=0;
PARAB=0;

    if INeedCircle
        params.grid.Xcenter=(params.grid.XLEFT+params.grid.XRIGHT)/2;
        params.grid.Ycenter=(params.grid.YLEFT+params.grid.YRIGHT)/2;
        circtol=1e-12;
        MhFermi=@(r,R,kT) 1./ ( 1+exp(( r-R )/kT )); %double(E-EF <1.5*kT)*
        R=params.grid.Lx*0.25; %radius of a cylinder

        Sth=pi*R^2; %theoretical area of the cylynder crossection
        kT=R*9e-3/(2^(params.grid.dx-5)); %smoothening parameter for the boundary
    %     kT=R*2e-3;
        [XX,YY]=meshgrid(params.grid.XLEFT+ params.grid.hx/2+...
            (0:(2^params.grid.dx-1))*params.grid.hx,params.grid.YLEFT+ params.grid.hy/2+(0:(2^params.grid.dy-1))*params.grid.hy);
        %___Geometry___
        Rs=sqrt(((XX-params.grid.Xcenter)).^2+((YY-params.grid.Ycenter)).^2);
        relerrS(i)=1;
        
        count=0;
        maxcount=500;
        while (relerrS(i)>circtol)&(count<maxcount)
            count=count+1;
            ffmask=MhFermi(Rs,R,kT);
            ffchi= (params.parts.eps-1)*params.rhs.k0^2*ffmask(:);
            S=sum(sum(ffmask))*params.grid.hx*params.grid.hy;
            oversize(i)=abs(S)/abs(Sth);
            relerrS(i)=abs(S-Sth)/abs(Sth);
            R=R*sqrt(1/oversize(i));
            kT=kT*sqrt(1/oversize(i));
        end
        
        
     
%(optional) Multi-scale material stuff    
% % % % % % % % % % % %         lcubex=0.4;
% % % % % % % % % % % %         lcubey=0.4;   
% % % % % % % % % % % %         params1=params;
% % % % % % % % % % % %         params1.dx=dx-log2(size(macromask,1));
% % % % % % % % % % % %         params1.dy=dy-log2(size(macromask,2));      
% % % % % % % % % % % %         Nx1=2^params1.dx;
% % % % % % % % % % % %         Ny1=2^params1.dy;
% % % % % % % % % % % %         Nz1=1;     
% % % % % % % % % % % % 
% % % % % % % % % % % %          ixmin=ceil(Nx1*(.5-lcubex))+1; %this is basically analytical mapping from 0..1 coordinate to index number
% % % % % % % % % % % %          ixmax=ceil(Nx1*(.5+lcubex));
% % % % % % % % % % % %          iymin=ceil(Ny1*(.5-lcubey))+1;
% % % % % % % % % % % %          iymax=ceil(Ny1*(.5+lcubey));
% % % % % % % % % % % %          material_mask= geom_box(params1,ixmin,ixmax,iymin,iymax,0,Nz1);
% % % % % % % % % % % %          
% % % % % % % % % % % %          
% % % % % % % % % % % %         material_mask=phize(material_mask,macromask,params1)       


        
    elseif PARAB %%parabolic mirror
        params.grid.Xcenter=(params.grid.XLEFT+params.grid.XRIGHT)/2;
        params.grid.Ycenter=(params.grid.YLEFT+params.grid.YRIGHT)/2;
        [XX,YY]=meshgrid(params.grid.XLEFT+ params.grid.hx/2+...
            (0:(2^params.grid.dx-1))*params.grid.hx,params.grid.YLEFT+ params.grid.hy/2+(0:(2^params.grid.dy-1))*params.grid.hy);
        smooth=@(bin)(1+tanh(-params.grid.Nx*bin/2))/2;
        R1=3;
        R2=R1;
        L1=sqrt(R1^2-1);
        L2=sqrt(R2^2-1);
        
        
         fun1=@(x,y)(x-L1-1+(R1-L1)).^2+y.^2-R1^2;
         fun2=@(x,y)(x+L2-1+(R2-L2)).^2+y.^2-R2^2;
        
%         fun=@(x,y)abs(y)-0.25
        
%         fun=YY+0.75-XX.^2;
%          ffmask=smooth(fun);
          ffmask=smooth(fun1(XX,YY)).*smooth(fun2(XX,YY));

        ffmask=ffmask.';
        figure;imagesc(ffmask);axis equal tight;colorbar
        
        ffchi= (params.parts.eps-1)*params.rhs.k0^2*ffmask(:);
       keyboard
%% Manual boxes placement
 %box1
 elseif 1 %cube
        params.surfaces.xmin = -400/800*abs(params.grid.XLEFT);
        params.surfaces.xmax = 400/800*abs(params.grid.XRIGHT);
        params.surfaces.ymin = -400/800*abs(params.grid.XLEFT);
        params.surfaces.ymax = 400/800*abs(params.grid.XRIGHT);
        params.surfaces.ixmin=ceil(params.grid.Nx*(params.surfaces.xmin-params.grid.XLEFT)/params.grid.Lx)+1; %xmin, xmax are in nanometers
        params.surfaces.ixmax=ceil(params.grid.Nx*(params.surfaces.xmax-params.grid.XLEFT)/params.grid.Lx);
        params.surfaces.iymin=ceil(params.grid.Ny*(params.surfaces.ymin-params.grid.YLEFT)/params.grid.Ly)+1;
        params.surfaces.iymax=ceil(params.grid.Ny*(params.surfaces.ymax-params.grid.YLEFT)/params.grid.Ly);
        material_mask= geom_box(params,params.surfaces.ixmin,params.surfaces.ixmax,params.surfaces.iymin,params.surfaces.iymax,0,params.grid.Nz);
%         material_mask=ttj2;
        chi = (params.parts.eps-1)*params.rhs.k0^2*material_mask
    elseif 0    
   
        ttX=round(params.grid.XLEFT+params.grid.hx/2+params.grid.hx*tt_x(2,params.grid.dx),1e-14);%volume-centric lattice is assumed
        ttY=round(params.grid.YLEFT+params.grid.hy/2+params.grid.hy*tt_x(2,params.grid.dy),1e-14);
        ttXX=mtkron(ttX,tt_ones(2,params.grid.dy)); %correct!
        ttYY=mtkron(tt_ones(2,params.grid.dx),ttY);
        fermi=@(E,E0, kT) 1./( 1 + exp((E-E0)/kT));
        r=@(xy)sqrt(   ( xy(:,1)/params.grid.XLEFT ).^2+( xy(:,2)/params.grid.XLEFT ).^2   )
        Rlens=params.grid.Lx/4;
        
        % Luneburg lens distributed material 

        smooth=5e-2;
%         
%         epsfun=@(xy)fermi(r(xy),Rlens,Rlens*smooth).* (2-r(xy)/Rlens) + (1-fermi(r(xy),Rlens,Rlens*smooth))*1;
%          epsfun=@(xy)fermi(r(xy),Rlens,Rlens*smooth).* (2*Rlens./(r(xy)+1e-1)-1) + (1-fermi(r(xy),Rlens,Rlens*smooth))*1;

         epsfun=@(xy) 1+(params.parts.eps-1)*exp(-r(xy).^2/Rlens)
%         epsfun=@(xy)(2*Rlens./(r(xy)+1e-1)-1)
        
         material_mask=round(amen_cross({ttXX,ttYY},@(x)epsfun(x),params.tol*1e-1,...
         'trunc_method','svd','kickrank',100, 'nswp',50,'max_err_jumps',6,'zrank',10),params.tol);
        chi = params.rhs.k0^2*material_mask;
%          figure;imagesc(reshape(abs(full(material_mask)),[2^params.grid.dx,2^params.grid.dy]));colorbar;colormap(viridis(2000));
        
        %Eaton lens distributed material    
        %eps(r)=(2R)/r-1
    end

    
% %box2
%         xmin = 200/800*2;
%         xmax = 300/800*2;
%         ymin = 0;
%         ymax = 200/800*2;
%         ixmin=ceil(Nx*(xmin-XLEFT)/Lx)+1; %xmin, xmax are in nanometers
%         ixmax=ceil(Nx*(xmax-XLEFT)/Lx);
%         iymin=ceil(Ny*(ymin-YLEFT)/Ly)+1;
%         iymax=ceil(Ny*(ymax-YLEFT)/Ly);
%         material_mask= material_mask+geom_box(params,ixmin,ixmax,iymin,iymax,0,Nz);
    %% (optional)FFTfication
        if FFTSolverOn && ~INeedCircle && ~PARAB
        ffmask=reshape(full(material_mask),[2^params.grid.dy,2^params.grid.dx]);
%         ffmask=permute(ffmask,[2 1]); %without permute something is wrong
% % % % % % % %         ffchi= (params.parts.eps-1)*params.rhs.k0^2*ffmask(:);   
        ffchi= full(chi);
        end
%         axis equal tight;
        colorbar
    
%% --(legacy, optional)Sacred reshape for 3D to work
% ffmask=reshape(full(\material_mask),[2^dy,2^dx]);
% data=real(reshape(full(material_mask),[2^dz,2^dy,2^dx]));
% data=permute(data,[ 3 2 1]);
%______
%% QTT chi forming block
if QTTSolverOn
    %encoding full mask into tensor (needed in case of circle or loaded picture)
% material_mask = tt_tensor(reshape(ffmask,2*ones(1,d)),tol);
% or making it from tt material_mask
end
%% --Old 3D plotting module (I dont need this in the near future) 
if letmePlot3d&&(params.grid.d==3)
    [fullmask1_3D, X, Y, Z]  =  F1_LtoF1_3D(full(material_mask), params);
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
%% Source forming block
fexpX=round(tt_sin_cos(params.grid.dx,params.grid.hx*params.rhs.k0x,1)+1j*...
    tt_sin_cos(params.grid.dx,params.grid.hx*params.rhs.k0x,0),params.tol)*exp(1j*(params.rhs.k0x*params.grid.hx/2));
fexpY=round(tt_sin_cos(params.grid.dy,params.grid.hy*params.rhs.k0y,1)+1j*...
    tt_sin_cos(params.grid.dy,params.grid.hy*params.rhs.k0y,0),params.tol)*exp(1j*(params.rhs.k0y*params.grid.hy/2));
SOURCEFIELD1=mtkron(fexpX,fexpY);
clear fexpX fexpY;
b=SOURCEFIELD1;
%% FFT Solving block
if FFTSolverOn
%% drawing the mask
if DrawMasks
        figure;imagesc((abs(ffmask)));colormap(gray(2048));
end
    %% solution for total field
    Axh=@(x)Ax(Gdeltas_fft,ffchi,params,x)
    ffb=full(b);
%     Poit-like source
%     ffb=1j/4 * besselh( 0,2,params.rhs.k0*sqrt((XX+0.5).^2+YY.^2) ).';
%     ffb=ffb(:);
    disp('gmres pure')
    [MEMDATA,~]=memory;
    maxmem=3*10^9;%0.9*MEMDATA.MemUsedMATLAB
    % 7Gbyte of Memory is accessable on this machine
    maxbasissize=floor(maxmem/2^params.grid.d/8)
    maxbasissize=min(maxbasissize,400);
    tic
    [x_fft,flag,relres,iter,resvec]=gmres(Axh,ffb,maxbasissize,params.tol,50);
    time.fft(i)=toc
    data_fft{i}=reshape(x_fft,[2^params.grid.dx,2^params.grid.dy]);
    data_fft{i}=permute(data_fft{i},[2 1]);
%% solution for scattered field    
%     disp('gmres pure scattered only')
%     ffb=full(b)-Axh(full(b));
%     tic
%     [x_fft,flag,relres,iter,resvec]=gmres(Axh,ffb,400,tol,200);
%     toc
%% plotting GMRES convergence
PlotGMRESConvergence=0;
if PlotGMRESConvergence
     figure;plot(log10(resvec/norm(ffb)),'-o');
     xlabel('iter');
     ylabel('log10 of relres');
end
%% zero phase matching to Mie for comparison    
%     data_fft_buf=data_fft{i}; %size(data_fft_buf,1)/2 +
%     data_fft{i}=data_fft{i}*...
%     exp(  -1j*(  angle(data_fft_buf(1,1))+k0x*hx/2+k0y*hy/2  )  ); %this is pixel-precision phase alignment with Mie results
%     clear data_fft_buf
%     data1=reshape(data1(:,:,ceil(2^(dz-1))),[2^dx 2^dy]).*(1-ffmask).*double(R>=EF);
%% plotting GMRES field
if FFTSolverOn && DrawFFTSolutions
    figure
%     figure('name','GMRES FFT solution');
    pltfun=@(x)abs(x);
    if exist('XX')
        imagesc(XX(1,:),YY(:,1).',pltfun(data_fft{i}));
        hold on
        contour(XX,YY,abs(ffmask).',1,'LineWidth',2,'Color','k')
        contour(XX,YY,abs(ffmask).',1,'LineWidth',1,'Color','w')%     theta=linspace(0,2*pi,180);%     plot(Xcenter+.25*Lx*cos(theta), Ycenter+.25*Ly*sin(theta),'--w','linewidth', 1.0)
        hold off
    else
        imagesc(pltfun(data_fft{i}));
        hold on
        contour(abs(ffmask).',1,'LineWidth',2,'Color','k')
        contour(abs(ffmask).',1,'LineWidth',1,'Color','w')%     theta=linspace(0,2*pi,180);%     plot(Xcenter+.25*Lx*cos(theta), Ycenter+.25*Ly*sin(theta),'--w','linewidth', 1.0)
        hold off
    end
    maxintens=max(max(abs(data_fft{i})))
    colormap(magma(2000))
    caxis([0 maxintens])
    axis equal tight;
    colorbar;
    xlabel('x');
    ylabel('y');
    title('fft') 
end
[M,I]=max(abs(x_fft))
xmax=mod(I, params.grid.Nx);
ymax=(I-mod(I,params.grid.Nx))/params.grid.Ny;
% YY(ymax,1);

% f_theor=R1/2/ sqrt(params.parts.eps)-1 
% f_measured=1-R1+sqrt(R1^2-1)-XX(1,xmax)

end
%% Calculating Mie
 if CALC_MIE
    tic
    plotDielectricCylinderTotalFieldUnderPlaneWave
    toc
    data_mie{i}= Ez_tot_norm;%*exp(-1j*k0x*hx); %ADD FREAKING PHASE GAUGE
 end
%% Draw Mie

if DrawMIESolutions && (LOAD_MIE||CALC_MIE)
    drawfun=@(x)abs(x);
    figure('name','Mie solution');
    if exist('XX')
        imagesc(XX(1,:),YY(:,1).',drawfun(data_mie{i}.'));
        hold on
        contour(XX,YY,ffmask.',1,'LineWidth',2,'Color','k')
        contour(XX,YY,ffmask.',1,'LineWidth',1,'Color','w')%     theta=linspace(0,2*pi,180);%     plot(Xcenter+.25*Lx*cos(theta), Ycenter+.25*Ly*sin(theta),'--w','linewidth', 1.0)
        hold off
    else
        imagesc(drawfun(data_mie{i}.'));
        hold on
        contour(ffmask.',1,'LineWidth',2,'Color','k')
        contour(ffmask.',1,'LineWidth',1,'Color','w')%     theta=linspace(0,2*pi,180);%     plot(Xcenter+.25*Lx*cos(theta), Ycenter+.25*Ly*sin(theta),'--w','linewidth', 1.0)
        hold off
    end
    colormap(magma(2000))
    caxis([0 2])
    axis equal tight;
    colorbar;
    xlabel('y');
    ylabel('x');
    title('mie')
end

%% QTT Amen-solve block
    if QTTSolverOn
        disp('forming tt_A')
        tic
         if params.grid.D==2
             G=round(tt_qtoepl(ttGXY,[params.grid.dx,params.grid.dy]),params.tol);
%              G2=abinitioG(params);
        elseif D==3
            G=round(tt_qtoepl(ttGXYZ,[params.grid.dx,params.grid.dy,params.grid.dz]),params.tol*1e-2);
         end
                tic
        A=tt_eye(2,params.grid.d)+round(G*diag(chi),params.tol) % to avoid neglecting small parts of G due to rounding errors we add I separately
        
        disp('amen pure')
        
        x_qtt=amen_solve2(A,b,params.tol,'resid_damp',1.0,'kickrank',5,'max_full_size',10000,'trunc_norm','residual');
        time.qtt(i)=toc
        data_qtt{i}=x_qtt;
%     bscat=round(b-A*b,tol);
%     
        zcord=params.grid.Nz
%% -Reconstruction of outer field from the fields inside scatterers
%     temp_extfield=-round(G*round(chi.*(x_qtt+b),tol),tol); %NO DIAG() here !! every column is multiplied on its own material constant
%     
%      data_reconstr=real(reshape(full(temp_extfield),[2^dz,2^dy,2^dx]));
%      data_reconstr=permute(data_reconstr,[ 3 2 1]);
%      data_reconstr=reshape(data_reconstr(:,:,zcord),[2^dx 2^dy]);
%     figure('name','reconstructed full solution');imagesc(data_reconstr);colorbar;axis equal tight;
%     colormap(jet(2048)); 

     %caxis([-2,2])
%     hold on
%     contour(data,1,'LineWidth',4,'Color','k')
%     contour(data,1,'LineWidth',2,'Color','w')
%     hold off
%% --Full Field Downscale section
%     dwn_qtt=@(x)downscale_2d_qtt(x,params);
%     data_qtt{dx-experiment.dx0+1}=x_qtt;
%     pic_qtt=real(reshape(full(x_qtt),[2^dz,2^dy,2^dx]));
%     pic_qtt=permute(pic_qtt,[ 3 2 1]);
%     pic_qtt=reshape(pic_qtt(:,:,zcord),[2^dx 2^dy]);
%% --QTT Field Downscale section
%     kostyl-style tests for dwn_qtt
%     fqtt=abs(full(dwn_qtt(dwn_qtt(dwn_qtt(dwn_qtt(x_qtt))))));
%     fqtt=reshape(fqtt,[sqrt(size(fqtt,1)),sqrt(size(fqtt,1))]);
%     figure; imagesc(fqtt );axis equal tight;colorbar;    colormap(jet(2048)); caxis([0,2.5])
%% QTT Drawing section

        if DrawQTTSolutions && params.grid.d <23
            pic_qtt=abs(reshape(full(x_qtt),[2^params.grid.dx,2^params.grid.dy]));
            pic_qtt=permute(pic_qtt,[2 1]);
            figure('name','Amen solution');
            ffmask=reshape(full(material_mask),...
                [params.grid.Ny,params.grid.Nx]);
            if exist('XX')
                imagesc(XX(1,:),YY(:,1).',pic_qtt);
                hold on
                contour(XX,YY,ffmask.',1,'LineWidth',2,'Color','k')
                contour(XX,YY,ffmask.',1,'LineWidth',1,'Color','w')
                hold off
            else
                imagesc(pic_qtt);   
                hold on
                contour(real(ffmask).',1,'LineWidth',2,'Color','k')
                contour(real(ffmask).',1,'LineWidth',1,'Color','w')
                hold off
            end
            colorbar;
            colormap(magma(2000))
%             caxis(2*[-0 1])
            axis equal tight
            ylabel('x');
            xlabel('y');
            title('qtt') 
        end
    end
maxr(i)=max(ttGXY.r)
end

%test_MIE
%% Internal convergence test, downscaling, second norm(sort of, used frobenius norm for 2d instead)
% % % % % % 
% % % % % % % ups=@(x)upscale_2d(x);
% % % % % % dwn=@(x)downscale_2d(x);
% % % % % % if FFTSolverOn
% % % % % %     for i=1:size(data_fft,2)-1
% % % % % %         relerr_int(i)=norm(    dwn((data_fft{i+1}))-(data_fft{i}),'fro')/norm((data_fft{i}),'fro');
% % % % % %     end
% % % % % %     
% % % % % %     if CALC_MIE || LOAD_MIE % it looks like mie series is for params.parts.eps=1.2, R=?0.25
% % % % % %         for i=1:size(data_fft,2)
% % % % % %         relerr_ext(i)=norm(abs(data_mie{i})-abs(data_fft{i}),'fro')/norm(abs(data_fft{i}),'fro');
% % % % % %         end    
% % % % % %         figure;plot(experiment.dx0:experiment.dxmax,log2(relerr_ext),'o');xlabel('d');ylabel('log2 relerr_{ext}');axis equal tight;
% % % % % %     end
% % % % % %     
% % % % % % %     ord2=@(n)   log2(...
% % % % % % %     2*norm( data_fft{n} - dwn(data_fft{n+1}) ,'fro') / norm(  data_fft{n+1}-dwn(data_fft{n+2}) ,'fro')...
% % % % % % %     )
% % % % % % 
% % % % % %     if size(data_fft,2)-1>1
% % % % % %         figure('name','fft internal convergence');plot((experiment.dx0:experiment.dx0+(size(data_fft,2)-2)) + 1,log2(relerr_int),'o-');
% % % % % %         xlabel('dx');
% % % % % %         ylabel('log2(err_{int})');
% % % % % %     end
% % % % % % 
% % % % % % end
% % % % % % 
% % % % % % 
% % % % % % if QTTSolverOn && (experiment.dxmax-experiment.dx0)>0   
% % % % % %     dwn_qtt=@(x)downscale_2d_qtt(x,params);
% % % % % %     for i=1:size(data_qtt,2)-1
% % % % % %     relerr_qtt_int(i)=norm(  ( ( (dwn_qtt(data_qtt{i+1}))-(data_qtt{i})  ))  )/norm(( (data_qtt{i})));
% % % % % %     end
% % % % % %     % relerr
% % % % % % 
% % % % % % %     ord2=@(n)   log2(...
% % % % % % %         norm( ( data_qtt{n} )-dwn_qtt((data_qtt{n+1})) ) / norm(  dwn_qtt(data_qtt{n+1})-dwn_qtt(dwn_qtt(data_qtt{n+2})))...
% % % % % % %         )
% % % % % %     
% % % % % % %     relerr_qtt_int=@(n) norm(downscale_2d(data_qtt{2})-data_qtt{1},'fro')/norm(data_qtt{1},'fro')
% % % % % % 
% % % % % %     
% % % % % %     % ord3=@(n)   log2(norm(abs(data1{n})-dwn(abs(data1{n+1}))) /norm(   (abs(data1{n+1}))-(dwn(abs(data1{n+2})))) )
% % % % % %     figure;title('QTT internal convergence');
% % % % % %     plot((experiment.dx0:experiment.dx0+(size(data_qtt,2)-2)) + 1,log2(relerr_qtt_int),'o-');
% % % % % %     xlabel('dx');
% % % % % %     ylabel('log2 relative error');
% % % % % % %     axis equal tight;
% % % % % %     
% % % % % % end