%DOESNT WORK!!
function field=mie_cyl_die(r,Radius,eps,th,mie_b,mie_c,k0)
    k1=k0;
    k2=k0*sqrt(eps);
    Nharmonics=(numel(mie_c)-1)/2;
    field=zeros(size(r));
    
%     .*exp(1j*n*th)
    for n=-Nharmonics:Nharmonics
%         keyboard
        addition=double(r>=Radius) .*(-1j)^n*mie_c(n+Nharmonics+1).*besselh(n,2,k1*r).*exp(1j*n*th)+...
                 double(r<Radius)  .*(-1j)^n*mie_b(n+Nharmonics+1).*besselj(n  ,k2*r).*exp(1j*n*th);
        addition(isnan(addition))=0;
        field=field+addition;
    end
    

end