function y=gh_xy_LCN(xy,patch,params)
    %LCN with one-pixel patch
    k0=params.rhs.k0;
    hx=params.grid.hx;
    okrest=1e-2*hx; %this measure of proximity implies only self-term being replaced by the patch number
    %this is here for the future possibility of patching a higher adjacency
    %of singularity
    
    deltar2=  xy(:,1).^2+xy(:,2).^2;
    y=zeros(size(xy,1),1);
    sum(deltar2<=okrest)
    patch
    y(deltar2<=okrest) = patch;
    y(deltar2>okrest) =1j/4* besselh( 0,2,k0*sqrt( deltar2(deltar2>okrest)) );

end