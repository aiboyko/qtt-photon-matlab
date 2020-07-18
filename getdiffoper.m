function Dx=getdiffoper(Nx,h)
%creates a 1st derivative operator of 2nd order of precision. N->N

h=1;
Dx=sparse(Nx,Nx);

Dx=diag(sparse(ones(Nx-1,1)),1)-diag(sparse(ones(Nx-1,1)),-1);
Dx=Dx/(2*h);
% 
Dx(1,:)=sparse(1,Nx);
Dx(1,1)=(-3/2)/h;
Dx(1,2)=2/h;
Dx(1,3)=(-1/2)/h;

Dx(Nx,:)=sparse(1,Nx);
Dx(Nx,Nx)=(3/2)/h;
Dx(Nx,Nx-1)=(-2)/h;
Dx(Nx,Nx-2)=(1/2)/h;

end