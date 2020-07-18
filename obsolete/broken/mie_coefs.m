function [mie_b,mie_c]=mie_coefs(k0,R,eps,Nharmonics)
%dielectric version - PEC only
k1=k0;
k2=k0*sqrt(eps);
dbesselj= @(n,z)1/2 * (besselj(-1 + n, z) - besselj(1 + n, z));
dbesselh2=@(n,z)1/2 * (besselh(-1+n,2, z) - besselh(1+n,2, z));

    for n=-Nharmonics:Nharmonics
        mie_b(n+Nharmonics+1)= (  besselj(n,k1*R)*dbesselh2(n,k1*R) - dbesselj(n,k1*R)*besselh(n,2,k1*R)  )...
            /...
                               (  besselj(n,k2*R)*dbesselh2(n,k1*R) - dbesselj(n,k2*R)*besselh(n,2,k1*R)  );
                           
        mie_c(n+Nharmonics+1)= (- besselj(n,k2*R)*dbesselj(n,k1*R)  + dbesselj(n,k2*R)*besselj(n,k1*R)  )...
            /...
                               (  besselj(n,k2*R)*dbesselh2(n,k1*R) - dbesselj(n,k2*R)*besselh(n,2,k1*R)  );
    end
end