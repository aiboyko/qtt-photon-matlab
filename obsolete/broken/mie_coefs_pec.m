function [mie_b,mie_c]=mie_coefs_pec(k0,R,Nharmonics)
%lazy version - PEC only
    for n=-Nharmonics:Nharmonics
        mie_b(n+Nharmonics+1)=0;
        mie_c(n+Nharmonics+1)= -(besselj(n,k0*R))/besselh(n,2,k0*R);
    end
end