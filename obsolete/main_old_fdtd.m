%By Alexey Boyko, Skolkovo Institute of Science and Technology,
%a.boyko@skoltech.ru
%alexey.boyko@skolkovotech.ru
%Moscow, Russia, 2015

clear all;
close all;

%% Setting up general constants -------------------------------
tol=1e-3;

%% Setting up the grid parameters------------------------------
dx=6;dy=6;dz=6;
d=dx+dy+dz;

%% Setting up physical envinronment and geometry and calculating in on qtt lattice-------
eps=3.00;
Lx=1;Ly=1;Lz=1;
hx=Lx/(2^dx); hy=Ly/(2^dy);hz=Lz/(2^dz);

h=min([hx hy hz]);
tau=.5*h;


% creating the box
% bawks=str_geom_box(0.1*Lx,.901*Lx,0.1*Ly,.901*Ly,.01*Lz,.99*Lz);
% 
% %creating Nx-by-Ny array of round holes
% str_holes={};
% Nx=1;Ny=1;
% R=.15;
% for i=1:Nx
%     for j=1:Ny
%         str_holes = cat(2,str_holes,{str_geom_cyl_xy(...
%             [Lx*i/(Nx+1),...
%             Ly*j/(Ny+1),...
%             R],0)});
%     end
% end
% 
% %combining previous geometries together
% general=str2func(strcat( '@(xyz)',  strjoin(str_holes,'&&'),'&&',bawks ));
% %general=str2func(strcat( '@(xyz)', bawks ));
% 
% %general=@(xyz) 1;
% % part1=@(xyz)bawks(xyz)&&general(xyz);
% % % part1=@(xyz)cf_and(xyz, cat(2,{bawks},holes));
% % % part2=@(xyz) part1(xyz-[0 0 2^(dz-2)*1.3 ]);
% % all_parts=@(xyz) part1(xyz)+ part2(xyz);
% 
% %calculating the geometry on the qtt grid via amen_cross
% %Amen.
% % mask_of_coords=amen_cross(2*ones(d, 1),@(bin) general(binarr2xyz(bin-1,...
% %     dx,dy,dz,hx,hy,hz)), 1e-14,'nswp',40);
% 
% % [X,Y,Z]=meshgrid(hy*(0:(2^dy-1)),hx*(0:(2^dx-1)),hz*(0:(2^dz-1)));

%--
x1=4;x2=18;
y1=10;y2=24;
z1=10;z2=24;
outer_box=tt_box(dx,dy,dz,tol,x1,x2,y1,y2,z1,z2);%round(mtkron(tth3z,tth3y,tth3x),tol);

xh1=7;xh2=11;
yh1=14;yh2=18;
zh1=1;zh2=32;
hole=round(1-tt_box(dx,dy,dz,tol,xh1,xh2,yh1,yh2,zh1,zh2),tol);

xh1=30;xh2=31;
yh1=2;yh2=30;
zh1=2;zh2=30;
flash=tt_box(dx,dy,dz,tol,xh1,xh2,yh1,yh2,zh1,zh2);

mask_of_coords = tt_zeros(2,d);%round(outer_box.*hole,tol);
%--

%% Assemblying Maxwell operator------------------------------
nil=ttm_zeros(2,dx+dy+dz);

eps_Ex=round(mask_of_coords*(eps-1)+1,tol);
eps_Ey=round(mask_of_coords*(eps-1)+1,tol);
eps_Ez=round(mask_of_coords*(eps-1)+1,tol);

eps_Ex_m=diag(eps_Ex);
eps_Ey_m=diag(eps_Ey);
eps_Ez_m=diag(eps_Ez);
eps_All=round(tt_blockwise({eps_Ex_m,      nil,       nil ;
                            nil,           eps_Ey_m,  nil;
                            nil,           nil,       eps_Ez_m     }),tol);

%warning! mtkron takes tensor product is reverse order
deltaEx=mtkron(tt_eye(2,dz),                tt_eye(2,dy),                   1/hx * fdtd_tt_DE(2,dx,tol));
deltaEy=mtkron(tt_eye(2,dz),                1/hy * fdtd_tt_DE(2,dy,tol),    tt_eye(2,dx));
deltaEz=mtkron(1/hz * fdtd_tt_DE(2,dz,tol), tt_eye(2,dy),                   tt_eye(2,dx));

%so right now X is the largest harmonic, Z is the highest-freq harmonic

deltaHx=mtkron(tt_eye(2,dz),                tt_eye(2,dy),               1/hx * fdtd_tt_DH(2,dx,tol));
deltaHy=mtkron(tt_eye(2,dz),                1/hy * fdtd_tt_DH(2,dy,tol),tt_eye(2,dx));
deltaHz=mtkron(1/hz * fdtd_tt_DH(2,dz,tol), tt_eye(2,dy),               tt_eye(2,dx));
% compressing deltas doesnt help much

curlE=tt_blockwise({nil,           -deltaEz,       deltaEy ;
                    deltaEz,        nil,            -deltaEx;
                    -deltaEy,       deltaEx,        nil     });

curlH=tt_blockwise({nil,           -deltaHz,       deltaHy ;
                    deltaHz,        nil,            -deltaHx;
                    -deltaHy,       deltaHx,        nil     });   

nil_nabla=tt_blockwise({nil,nil,nil;
                        nil,nil,nil;
                        nil,nil,nil}); 

                    
% AE2H=round(-tau*nablaE*eps_All,tol);
AE2H=round(-tau*curlE,tol);
AH2E=round(tau*curlH,tol);

% max(abs(eig(full(AE2H))))
% max(abs(eig(full(AH2E))))
                    
% A=tt_blockwise({nil_nabla,  -nablaE*eps_All;
%                 nablaH,  nil_nabla});
%             mem(A)
% A=round(A,tol)     ;
% mem(A)

%% converting geometry into FULL format and visualising -------------

[fullmask_reshaped,X,Y,Z]=QTTFIELDto3DFIELD(mask_of_coords,dx,dy,dz,hx,hy,hz);

letmePlot3d=1;
if letmePlot3d

figure(1)
iso_mask=isosurface(X,Y,Z,fullmask_reshaped);
p=patch(iso_mask);
%p.FaceAlpha=0.6;
p.EdgeColor = 'none';
daspect([1,1,1])
view(3); axis tight
p.FaceColor = 'green';
camlight 
lighting gouraud

keyboard
end
%% Time integration-----------------------------

maxt=1000;
%ICS for t=0
%at t=0
icsHx=5*flash;
icsHy=tt_zeros(2,d);
icsHz=tt_zeros(2,d);

icsH=tt_vectcat_ext(icsHx,icsHy,icsHz);
icsH=round(icsH,tol);
HFIELD=icsH;
movie_HFIELD{1}=HFIELD;

%t=1/2 tau
icsEx=tt_zeros(2,d);
icsEy=tt_zeros(2,d);
icsEz=5*flash;%tt_ones(2,d);
icsE=tt_vectcat_ext(icsEx,icsEy,icsEz);
icsE=round(icsE,tol);
EFIELD=icsE;

movie_EFIELD{1}=EFIELD;
Eframe=tt_vectsplit_ext(EFIELD,3);
Eframe3d=QTTFIELDto3DFIELD(Eframe,dx,dy,dz,hx,hy,hz);
Eframe3dslice1(:,:,1)=Eframe3d(:,:,2^(dz-1));
% Eframe3dslice2(:,:,1)=reshape(Eframe3d(:,15,:),[2^dx,2^dz]);


%icsFIELD=tt_vectcat_ext(icsE,icsH);

%Yee cell implies 2 time-staggered grids for E (fractional timesteps) and B
%(integer timesteps)
%because of seeming difficulty of updating only half of the ttvect, we
%shall consider a block-wise(separate for E and B) time-updating procedure

fig1=figure(2);
ax1=axes();
axis equal tight
set(ax1,'nextplot',...
'replacechildren');

% fig2=figure(3);
% ax2=axes();
% axis equal tight
% set(ax2,'nextplot',...
% 'replacechildren');

fullmask_reshaped_slice1=fullmask_reshaped(:,:,2^(dz-1));
fullmask_reshaped_slice2=reshape(fullmask_reshaped(:,15,:),[2^dy,2^dz]);

for t=1:maxt
    t
    HFIELD=HFIELD+mvk3(AH2E,HFIELD,tol);
    HFIELD=round(HFIELD,tol);
    movie_HFIELD{t+1}=HFIELD;
    
    EFIELD=EFIELD+mvk3(AE2H,EFIELD,tol);
    EFIELD=round(EFIELD,tol);
    movie_EFIELD{t+1}=EFIELD;
    
    %experiment-show on the fly
    Eframe=tt_vectsplit_ext(EFIELD,3);
    Eframe3d=QTTFIELDto3DFIELD(Eframe,dx,dy,dz,hx,hy,hz);
    
    Eframe3dslice1(:,:,t+1)=Eframe3d(:,:,2^(dz-1));
%    
    Eframe3dslice2(:,:,t+1)=reshape(Eframe3d(:,15,:),[2^dy,2^dz]);
  
    figure(2)
%     contourf(Eframe3dslice1(:,:,t),'Edgecolor','none');
    pcolor(Eframe3dslice1(:,:,t));
    shading interp
    hold on
    contour(fullmask_reshaped_slice1 ,[2],'LineWidth',1,'LineColor','black');
    axis equal tight
    hold off
    
    figure(3)
%     contourf(Eframe3dslice2(:,:,t),'Edgecolor','none');
    pcolor(Eframe3dslice2(:,:,t));
    shading interp
    hold on
    contour(fullmask_reshaped_slice2,[2],'LineWidth',1,'LineColor','black');
    axis equal tight
    hold off
    
    drawnow;
end

%--------------------------------------------------------------------------
%% Data Movie-Making----------------------------
%--------------------------------------------------------------------------
% for tt=1:t
%     Eframe=movie_EFIELD{tt};
%     Eframe=tt_vectsplit_ext(Eframe,3);
%     Eframe3d=QTTFIELDto3DFIELD(Eframe,dx,dy,dz,hx,hy,hz);
%     Eframe3dslice1(:,:,tt)=Eframe3d(:,:,2^(dz-1));
% end

%--------------------------------------------------------------------------
%-----------------------------Field Animation-----------------------------
%--------------------------------------------------------------------------
% framestep=1;
% fig1=figure(2);
% ax1=axes();
% set(ax1,'nextplot',...
% 'replacechildren');
% % obj=contourf(Eyframe3dslice,'Edgecolor','none')
% % deltaXYZ_v=deltaXYZ('Ez').*[hx hy hz];
% fig2=figure(3);
% ax2=axes();
% set(ax2,'nextplot',...
% 'replacechildren');
% 
% for tt=1:framestep:t
%     figure(fig1)
%     contourf(X(:,:,2^(dz-1)),Y(:,:,2^(dz-1)),Eframe3dslice1(:,:,tt),'Edgecolor','none');
%     hold on
%     contour(X(:,:,2^(dz-1)),Y(:,:,2^(dz-1)),fullmask_reshaped(:,:,16),[2],'LineWidth',1,'LineColor','black')
%     axis equal tight
%     hold off
%         
%     pause(0.005);
%     drawnow;
% end

 V=full(eps_Ex);
 HA=vol3d('cdata',V,'texture','3D');
 view(3);  
 axis tight;  daspect([1 1 1])

