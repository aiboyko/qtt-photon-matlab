function y=G_full_helm3d(params)
    
    k0=params.rhs.k0;
    
    Nx = 10;%params.grid.Nx;
    Ny = Nx;%params.grid.Ny;
    Nz = Nx;%params.grid.Nz;
    Lx = params.grid.Lx;
    Ly = params.grid.Ly;
    Lz = params.grid.Lz;
    R = sqrt(Lx^2 + Ly^2 + Lz^2)+0.01;
    
    overs=4; %=4 in the paper
    ksx = 2*pi / (Lx*overs) * (-overs*floor(Nx/2):overs*floor(Nx/2) -1 );
    ksy = 2*pi / (Ly*overs) * (-overs*floor(Ny/2):overs*floor(Ny/2) -1 );
    ksz = 2*pi / (Lz*overs) * (-overs*floor(Nz/2):overs*floor(Nz/2) -1 );

    [kksx,kksy,kksz] = meshgrid(ksx,ksy,ksz);

    K_vals = helm3d_kernel(sqrt(kksx.^2 + kksy.^2 + kksz.^2),k0,R);

    [i1,i2,i3] = meshgrid(0:size(K_vals,1)-1,0:size(K_vals,2)-1,0:size(K_vals,3)-1);

    y = ifftn(K_vals).*(-1).^(i1+i2+i3);


    y = y(1:Nx,1:Ny,1:Nz);

end