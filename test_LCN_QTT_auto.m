%from mathematica
%  -8.4656627077126*10^-7 + 2.3841708304603*10^-7 I
tol=5e-14;
for i=22:22
ddx=i;
ddy=ddx;
NNx=2^ddx;
NNy=2^ddy;
hhx = params.grid.hx / NNx;
hhy = params.grid.hy / NNy;
micro_Xdelta=round(-params.grid.hx/2 + hhx/2 + hhx*tt_x(2,ddx),1e-14);
micro_Ydelta=round(-params.grid.hy/2 + hhy/2 + hhy*tt_x(2,ddy),1e-14);


ttmicroXX=mtkron(micro_Xdelta,tt_ones(2,ddy));
ttmicroYY=mtkron(tt_ones(2,ddx),micro_Ydelta);

ghh=@(xy)hhx*hhy*1j/4*...
    besselh( 0,2,params.rhs.k0 * sqrt(  xy(:,1).^2 + xy(:,2).^2  ));

quad_field=round(amen_cross({ttmicroXX,ttmicroYY},@(xy)ghh(xy),tol*1e-1,...
        'trunc_method','svd','kickrank',100, 'nswp',100,'max_err_jumps',15,'zrank',17),tol);
% figure; imagesc(real(reshape( full(quad_field),[NNx,NNy] )))
patch(i)=sum(quad_field)
relpatch(i)=abs((patch(i)-patch(1))/patch(1));
end