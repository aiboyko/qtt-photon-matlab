function y=gh_xxyy_Nyst_selfheal(xxyy,patch,params)
    %tensor for full matrix (not entity matrix) LCN with one-pixel patch
    k0=params.k0;
    hx=params.hx;
    hy=params.hy;
    
    okrest=1e-3*hx;
    deltar2=(  xxyy(:,1)-xxyy(:,2)  ).^2+(  xxyy(:,3)-xxyy(:,4)  ).^2;
    y=zeros(size(xxyy,1),1);
    y(deltar2<=okrest) = hx*hy*1j/4* besselh( 0,2,k0*sqrt(hx*hy)*params.quadratures.softeningcoef );
    y(deltar2>okrest)  = hx*hy*1j/4* besselh( 0,1,k0*sqrt( deltar2(deltar2>okrest)) );
end