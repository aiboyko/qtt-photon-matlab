function y=gh_xy_Nyst_shifted(xy,params)
    %LCN with one-pixel patch
    k0=params.rhs.k0;
    hx=params.grid.hx;
    hy=params.grid.hy;
       
    deltar2=  (  xy(:,1)+params.grid.hx/2  ).^2+xy(:,2).^2;
    y=zeros(size(xy,1),1);
    y = 1j/4* besselh( 0,2,k0*sqrt( deltar2) );

end