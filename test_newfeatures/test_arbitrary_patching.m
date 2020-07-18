%macroharmonic (cores that are responsible for large jumps) is far right
% 

angle=0;
zoom=1;

I = (255-rgb2gray(imread('jackel2.png')));
I=imrotate(I,angle,'bilinear','crop');
I=imresize(I,[128 128]);
tol=1e-5;
SIZ=size(I,1);
Z=zeros(SIZ);
Z(1:size(I,1),1:size(I,2))=double((I(:,:,1)));

%Z=ones(4)
patch=reshape(Z,2*ones(1,log2(numel(Z))));
ttpatch=tt_tensor(patch)
dx_patch=log2(numel(patch))/2;
dy_patch=log2(numel(patch))/2;
d_patch=dy_patch+dx_patch
% 

% lito=zeros(8);
% lito(7,1)=1;
% lito(1,3)=1;figure;
% dx_lito=log2(numel(lito))/2;
% dy_patch=log2(numel(patch))/2;
% 
% ttpiss=tt_tensor(reshape(lito,2*ones(1,log2(numel(lito)))),1e-5)
%  
%  % compp=kron(piss, reshape(patch,[ 4 4]))
% compt=tkron( ttp,ttpiss);
% compt=permute(compt,[ 1 2 5 6 7 3 4 8 9 10],...
%     1e-5);
% patch=reshape(ones(4),2*ones(1,4));
% dx_patch=log2(numel(patch))/2;
% dy_patch=log2(numel(patch))/2;
% d_patch=dy_patch+dx_patch
% 
% ttp=tt_tensor(patch)



dx_lito=6;
dy_lito=6;
lito=zeros(2^dx_lito,2^dy_lito);
lito(1,1)=1;

d_lito=dy_lito+dx_lito;

dx=dx_lito+dx_patch;
dy=dy_lito+dy_patch;
d=dx+dy;

params0.dx=dx;
params0.dy=dy;
params0.d=d;
params0.tol=tol;

ttlito=tt_tensor(reshape(lito,2*ones(1,d_lito)),tol);
 
 % compp=kron(piss, reshape(patch,[ 4 4]))
compt=tkron( ttpatch,ttlito);
compt=permute(compt,[ 1:dx_patch d_patch+1:d_patch+dx_lito dx_patch+1:d_patch  d_patch+dx_lito+1:d_patch+d_lito],...
    1e-5);


s=0;
for i=0:(2^dy_lito)-1
    
% xshift=(2^(dx-1)-2^(dx_patch-1));
% yshift=(2^(dy-1)-2^(dy_patch-1));
xshift=0;
yshift=min((2^dy_patch)*i,2^dy-zoom*2^dy_patch)+10*rand(1);

cx1=ceil(xshift)-xshift;
cx2=1-cx1;
cy1=ceil(yshift)-yshift;
cy2=1-cy1;

s=round(s+qtt_2d_shift(compt,xshift,yshift,params0),tol);
end
s

fc=reshape(full(s),[ 2^dx 2^dy ]);
% 
% figure(1);imagesc(fc);colorbar;axis equal
