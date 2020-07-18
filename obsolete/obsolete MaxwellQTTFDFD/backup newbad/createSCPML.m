function [invsx,invsy,invsz]=createSCPML(params,pml)

pmltol=1e-8;
    dx=params.dx;
    dy=params.dy;
    dz=params.dz;
    
    hx=params.hx;
    hy=params.hy;
    hz=params.hz;
    
    Lx=params.Lx;
    Ly=params.Ly;
    Lz=params.Lz;
    tol=params.tol;
    k0=params.k0;
    omega=params.omega;
    eps0=params.eps0;
    mu0=params.mu0;
    etha0=params.etha0;
    

    xlowwidth=pml.xlowwidth;
    xhighwidth=pml.xhighwidth;
    ylowwidth=pml.ylowwidth; 
    yhighwidth=pml.yhighwidth;
    zlowwidth=pml.zlowwidth;   
    zhighwidth=pml.zhighwidth;
    
    d=dx+dy+dz;
    
%     % 1e0 -+
%1e6 reflections
%1e9 refl
    R=pml.R;
    m=pml.m;
    sigma.xhigh=-(m+1)*log(R)/(2*etha0*xhighwidth);
    sigma.xlow=-(m+1)*log(R)/(2*etha0*xlowwidth);
    sigma.yhigh=-(m+1)*log(R)/(2*etha0*yhighwidth);
    sigma.ylow=-(m+1)*log(R)/(2*etha0*ylowwidth);
    sigma.zlow=-(m+1)*log(R)/(2*etha0*zlowwidth);
    sigma.zhigh=-(m+1)*log(R)/(2*etha0*zhighwidth);
    
    pml_xhigh=Lx - xhighwidth;
    pml_xlow = xlowwidth;
    pml_yhigh=Ly - yhighwidth;
    pml_ylow = ylowwidth;
    pml_zhigh=Lz - zhighwidth;
    pml_zlow = zlowwidth;
    
     sxf=@(xyz) ...
              (xyz(1) >=  pml_xhigh).*(1 + 1j* abs(( xyz(1)-pml_xhigh)/xhighwidth  ).^m*sigma.xhigh/(eps0*omega))+...
               (xyz(1) <=  pml_xlow ).*(1 + 1j* abs(( xyz(1)-pml_xlow )/xlowwidth   ).^m*sigma.xlow /(eps0*omega))+...
               (xyz(1) > pml_xlow ).*( xyz(1) <pml_xhigh );
 
 
%                 (xyz(1) >  pml_xhigh).*(   1 - 1j.*3.2e+1* abs(( xyz(1)-pml_xhigh)/xhighwidth  ).^m  )+...
%                (xyz(1) <  pml_xlow ).*(1 - 1j.*3.2e+1* abs(( xyz(1)-pml_xlow )/xlowwidth   ).^m    )+...
%                (xyz(1) >= pml_xlow ).*( xyz(1) <=pml_xhigh );



kt=sigma.xhigh/(eps0*omega)

%                 (xyz(1)>pml_xhigh)*(1+1j*sigmamax/k0*  sin(pi/2*((xyz(1)-pml_xhigh)/xhighwidth))^2)*(1 + (amax-1)*((xyz(1) - pml_xhigh)/ xhighwidth)^p ) +...
%                 (xyz(1)<pml_xlow )*(1+1j*sigmamax/k0*  sin(pi/2*((pml_xlow -xyz(1))/xlowwidth))^2) *(1 + (amax-1)*((xlowwidth - xyz(1))/ xlowwidth )^p ) +...
%                 (xyz(1)>=pml_xlow)*(xyz(1)<=pml_xhigh) ;s
...

    syf=@(xyz) 1;        
%                (xyz(2) >= pml_yhigh).*(1+1j*a*abs((xyz(2)-pml_yhigh)/yhighwidth).^p)+...
%                (xyz(2) <= pml_ylow ).*(1+1j*a*abs((xyz(2)-pml_ylow )/ylowwidth ).^p)+...
%                (xyz(2) >  pml_ylow ).*(xyz(2)<pml_yhigh);
% %     (xyz(2)>pml_yhigh)*(1+1j*sigmamax*  sin(pi/2*((xyz(2)-pml_yhigh)/ylowwidth))^2)*(1 + (amax-1)*((xyz(2)  -  pml_yhigh)/yhighwidth)^p ) +...
%                 (xyz(2)<pml_ylow)*(1-1j*sigmamax*  sin(pi/2*((pml_ylow -xyz(2))/ylowwidth))^2)*(1 + (amax-1)*((ylowwidth-xyz(2))/ ylowwidth)^p ) +...
%                (xyz(2)>=pml_ylow)*(xyz(2)<=pml_yhigh) ;
    szf=@(xyz) 1;
%     (xyz(3)>pml_zhigh)*(1-1j*sigmamax*  sin(pi*((xyz(3)-pml_zhigh)/(2*zlowwidth))^2))*(1 + (amax-1)*((xyz(3)  -  pml_zhigh)/zhighwidth)^p ) +...
%                 (xyz(3)<pml_zlow)*(1-1j*sigmamax*  sin(pi*((pml_zlow -xyz(3))/(2*zlowwidth))^2))*(1 + (amax-1)*((zlowwidth-xyz(3))/ zlowwidth)^p ) +...
%                (xyz(3)>=pml_zlow)*(xyz(1)<=pml_zhigh) ;
           
    invsxf=@(xyz) 1/sxf(xyz);
    invsyf=@(xyz) 1/syf(xyz);
    invszf=@(xyz) 1/szf(xyz);
    
invsx=  round(amen_cross(2*ones(d, 1),@(bin) invsxf(bin2xyz(bin-1,...
     params)), pmltol,'nswp',40),tol);  
invsy=  round(amen_cross(2*ones(d, 1),@(bin) invsyf(bin2xyz(bin-1,...
     params)), pmltol,'nswp',40),tol); 
invsz=  round(amen_cross(2*ones(d, 1),@(bin) invszf(bin2xyz(bin-1,...
     params)), pmltol,'nswp',40),tol); 
 
   
       invsx_full=F1_LtoF1_3D(full(invsx), params);
       figure();imagesc(real( invsx_full(:,:,1)));colorbar; 
       figure();imagesc(imag( invsx_full(:,:,1)));colorbar;
figure;plot(real( invsx_full(:,:,1)))
       figure;plot(imag( invsx_full(:,:,1)))
keyboard
      
end