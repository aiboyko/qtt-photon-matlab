
params.grid.XLEFT=-.5; %use -2lambda, +2 lambda for comparison with Mie
params.grid.XRIGHT=.5;
params.grid.YLEFT=-.5;
params.grid.YRIGHT=.5;

params.grid.Lx=params.grid.XRIGHT-params.grid.XLEFT; 
params.grid.Ly=params.grid.YRIGHT-params.grid.YLEFT;
params.grid.Lz=1;

% |[ . ][ . ][ . ][ . ]|
% |<--------L--------->|

params.grid.Nx=11;
params.grid.Ny=params.grid.Nx;
