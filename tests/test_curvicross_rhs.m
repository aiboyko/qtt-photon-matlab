
xleft=-1;
xright=1;
yleft=-1;
yright=1;

Rmax=1;


Lx=xright-xleft;
Ly=yright-yleft;
dx=5;
dy=dx;
d=dx+dy;

Nx=2^dx;
Ny=2^dy;

hx=Lx/Nx;
hy=Ly/Ny;
ttX=hx*(tt_x(2,dx+1))+hx/2;%volume-centric lattice is assumed [ . ][ . ]
ttY=hy*(tt_x(2,dy+1))+hy/2;

dr=dx;
dphi=dx;
NR=2^dr;
NPHI=2^dphi;

hr=Rmax/NR;
hphi=2*pi/NPHI;

ttR=hx*(tt_x(2,dx+1))+hx/2;%volume-centric lattice is assumed [ . ][ . ]
ttPHI=2*hphi*(tt_x(2,dy+1))+hy/2;

mej=params.quadratures.maxerrorjumps;
k0x=params.rhs.k0x;
k0y=params.rhs.k0y;
tol=params.tol;
softeningcoef=params.quadratures.softeningcoef;

%general approach with composition of functions
%gh_curv2xy=@(coord) [coord(:,1).*cos(coord(:,3)),coord(:,2).*cos(coord(:,4)), coord(:,1).*sin(coord(:,3)) , coord(:,2).*sin(coord(:,4))];

%direct computation with cosine theorem
%gh_total=@(coord) 1./sqrt( coord(:,1).^2+coord(:,2).^2 - 2*coord(:,1).*coord(:,2).*cos(coord(:,3)-coord(:,4))  +   tol*1e1)

%generating standard 2D mesh
ttXX=mtkron(ttX,tt_ones(2,dy)); %correct!
ttYY=mtkron(tt_ones(2,dx),ttY);

ttRR=mtkron(ttR,tt_ones(2,dy)); %correct!
ttPHIPHI=mtkron(tt_ones(2,dx),ttPHI);


% gh=@(xy)   xy(:,1).^2+xy(:,2).^2  ;

rphi2xy=@(rphi)[rphi(:,1).*cos(rphi(:,2)) , rphi(:,1).*sin(rphi(:,2))];
gh=@(xy) exp(1j*( k0x*xy(:,1)+k0y*xy(:,2) ));



ttfun=round(amen_cross({ttXmesh,ttYY},@(x)gh(x),tol*1e-1,...
'trunc_method','svd','kickrank',100, 'nswp',100,'max_err_jumps',9,'zrank',13),tol);

