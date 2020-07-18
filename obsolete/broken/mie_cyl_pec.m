function field=mie_cyl_pec(r,Radius,th,mie_b,mie_c,k0)
    Nharmonics=(numel(mie_c)-1)/2;
    field=zeros(size(r));
    
    for n=-Nharmonics:Nharmonics
%         keyboard
        addition=(-1j)^n*mie_c(n+Nharmonics+1).*double(r>Radius).*besselh(n,2,k0*r).*exp(1j*n*th);
        addition(isnan(addition))=0;
        field=field+addition;

    end
    

end