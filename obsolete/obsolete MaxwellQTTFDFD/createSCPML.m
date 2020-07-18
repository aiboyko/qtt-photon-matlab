function [invsx,invsy,invsz]=createSCPML(params,pml)

    dx=params.dx;
    dy=params.dy;
    dz=params.dz;
    hx=params.hx;
    hy=params.hy;
    hz=params.hz;
    tol=params.tol;
    k0=params.k0;

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
    a.xhigh=-(m+1)*log(R)/(2*xhighwidth*k0)
    a.xlow=-(m+1)*log(R)/(2*xlowwidth*k0);
    a.yhigh=-(m+1)*log(R)/(2*yhighwidth*k0);
    a.ylow=-(m+1)*log(R)/(2*ylowwidth*k0);
    a.zlow=-(m+1)*log(R)/(2*zlowwidth*k0);
    a.zhigh=-(m+1)*log(R)/(2*zhighwidth*k0);
    
    pml_xhigh=(2^dx)-1 - xhighwidth;
    pml_xlow = xlowwidth;
    pml_yhigh=(2^dy)-1 - yhighwidth;
    pml_ylow = ylowwidth;
    pml_zhigh=(2^dz)-1 - zhighwidth;
    pml_zlow = zlowwidth;

    sxf=@(xyz) ...
               (xyz(1) >= pml_xhigh).*(1+1j*a.xhigh*abs((xyz(1)-pml_xhigh)/xhighwidth).^m)+...
               (xyz(1) <= pml_xlow ).*(1+1j*a.xlow*abs((xyz(1)-pml_xlow)/xlowwidth ).^m)+...
               (xyz(1) >  pml_xlow ).*(xyz(1)<pml_xhigh);

        
%                         (xyz(1)>pml_xhigh)*(1+1j*sigmamax/k0*  sin(pi/2*((xyz(1)-pml_xhigh)/xhighwidth))^2)*(1 + (amax-1)*((xyz(1) - pml_xhigh)/ xhighwidth)^p ) +...
%                 (xyz(1)<pml_xlow )*(1-1j*sigmamax/k0*  sin(pi/2*((pml_xlow -xyz(1))/xlowwidth))^2) *(1 + (amax-1)*((xlowwidth - xyz(1))/ xlowwidth )^p ) +...
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
     params)), 1e-10,'nswp',40),tol);  
invsy=  round(amen_cross(2*ones(d, 1),@(bin) invsyf(bin2xyz(bin-1,...
     params)), 1e-10,'nswp',40),tol); 
invsz=  round(amen_cross(2*ones(d, 1),@(bin) invszf(bin2xyz(bin-1,...
     params)), 1e-10,'nswp',40),tol); 
 
%     figure();
%     subplot(1,2,1);
%     plot(real(full(syz_x)))
%     subplot(1,2,2);
%     plot(imag(full(syz_x)))
% %     
%      invsx_full=F1_LtoF1_3D(full(invsx), dx,dy,dz,hx,hy,hz);
%      figure();imagesc(real( invsx_full(:,:,1)));colorbar; 
%      figure();imagesc(imag( invsx_full(:,:,1)));colorbar;
%      keyboard
end