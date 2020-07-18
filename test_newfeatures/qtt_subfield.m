function [y,newdx,newdy,newdz]=qtt_subfield (ttv,params,dxmax,dymax,dzmax)
%sparsifies the 3d-qtt tt_vector in dimensions X and Y
% freezes Z dimension in coordinate iz=1
dx=params.dx; 
dy=params.dy;
dz=params.dz;
d=dx+dy+dz;
newdx=min(dx,dxmax);
deltadx=dx-newdx;

newdy=min(dy,dymax);
deltady=dy-newdy;

newdz=min(dz,dzmax);
deltadz=dz-newdz;

if deltady>=1
    deltayvals=[ones(1,deltady-1) 2];
else
    deltayvals=[ ones(1,deltady-1)];
end
if deltadx>=1
    deltaxvals=[ones(1,deltadx-1) 2];
else
    deltaxvals=[ ones(1,deltadx-1)];
end

idxlist=[1: deltadz,dz+1:dz+deltady,dz+dy+1:dz+dy+deltadx];
idxvals=[2*ones(1,dz)  deltayvals deltaxvals   ];
%as one can see, the z value is put to 1 1 1 1 1 = 1
y=tt_subtensor(ttv,idxlist,idxvals);
end