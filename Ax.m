function mv=Ax(GA,chi,params,x)
    %2D Helmholtz matvec with (2N-1)*(2M-1) array of data, which is exactly what is needed to
    % assemble BTTB matrix with MxM block-Toeplitz  matrix with NxN-sized Toeplitz blocks
    %GA is a BTTB_entity in FFT format
    
    hx=params.grid.hx;
    hy=params.grid.hy;
    hz=params.grid.hz;

    dz=params.grid.dz;

    Nx=(size(GA,1)+1)/2;
    Ny=(size(GA,2)+1)/2;
    Nz=(size(GA,3)+1)/2;

    x1=chi.*x;
    if dz==0
        z=zeros(2*Nx-1,2*Ny-1);
        z(1:Nx,1:Ny)=reshape(x1,[Nx Ny]);
        mv=ifftn(fftn(GA).*fftn(z));
        mv=reshape(mv(1:Nx,1:Ny),[Nx*Ny,1]);
        mv=x+mv;
    else
        z=zeros(2*Nx-1,2*Ny-1,2*Nz-1);
        z(1:Nx,1:Ny,1:Nz)=reshape(x1,[Nx Ny Nz]);
        mv=ifftn(fftn(GA).*fftn(z));
        mv=reshape(mv(1:Nx,1:Ny,1:Nz),[Nx*Ny*Nz,1]);
        mv=x+mv;
    end
end