function vargout=F3L_visualize(F3L,params,varargin)
%Draws a field from linear non-tensor data

field='Ez';

dx=params.dx;
dy=params.dy;
dz=params.dz;
hx=params.hx;
hy=params.hy;

N=2^(dx+dy+dz);
%pre-coded zcoord of out vis
z_coord=2^(dz-1);
% fun='re';

for k = 1:length(varargin)
    if strcmp(varargin{k},'field')
        field = varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    elseif strcmp(varargin{k},'z_coord')
        z_coord = varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    elseif strcmp(varargin{k},'fig')
        fig = varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    elseif strcmp(varargin{k},'func')
        func = varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];    
    elseif strcmp(varargin{k},'name')
        name = varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];          
    end
end

if field == 'Ex'
    indexrange=(1):(N);
elseif field == 'Ey'
    indexrange=(N+1):(2*N);
elseif field == 'Ez'
    indexrange=(2*N+1):(3*N);
elseif field == 'Hx'
    indexrange=(3*N+1):(4*N);
elseif field == 'Hy'
    indexrange=(4*N+1):(5*N);
elseif field == 'Hz'
    indexrange=(5*N+1):(6*N);   
end

if ~exist('func')
    func='re'
end

if strcmp(func,'re')
    funh=@(x)real(x);
elseif strcmp(func,'im')
    funh=@(x)imag(x);
elseif strcmp(func,'abs')
    funh=@(x)abs(x);
elseif strcmp(func, 'arg')
    funh=@(x)angle(x); 
end

Fz_s_3d = F1_LtoF1_3D(F3L(indexrange),params);
if exist('fig')
fig=figure(fig);
else
fig=figure();
end

if exist('name')
fig.Name=name;
end

x=hx*(1:(2^dx-1));
y=hy*(1:(2^dy-1));

imagesc(x,y,funh(Fz_s_3d(:,:,2^(dz))));
title(strcat(func,'(',field,')'));

% axis equal tight;
colorbar;
if exist('params.mask3d')
hold on
contour(params.mask3d(:,:,z_coord),1,'color','Black');
hold off
end
vargout=fig;
end