function [shifted_ttv2d,shift_ttmatrix]=qtt_2d_shift(ttv2d,xshift,yshift,params)
dx=params.dx;
dy=params.dy;
tol=params.tol;

cx1=ceil(xshift)-xshift;
cx2=1-cx1;
cy1=ceil(yshift)-yshift;
cy2=1-cy1;

shift_ttmatrix=tkron(cx1*tt_shift(2,dx,floor(xshift))+...
cx2*tt_shift(2,dx,ceil(xshift)),...
    cy1*tt_shift(2,dy,floor(yshift))+...
    cy2*tt_shift(2,dy,ceil(yshift)) ...
    );
shifted_ttv2d=round(shift_ttmatrix*ttv2d,tol);
end