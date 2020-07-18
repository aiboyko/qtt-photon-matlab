function [syz_x,sxz_y,sxy_z,invsxz_y,invsyz_x,invsxy_z]=createUPML(params,pml)

%     dx=params.dx;
%     dy=params.dy;
%     dz=params.dz;
%     
%     hx=params.hx;
%     hy=params.hy;
%     hz=params.hz;
%     
%     Lx=params.Lx;
%     Ly=params.Ly;
%     Lz=params.Lz;
%     tol=params.tol;
%     k0=params.k0;
%     omega=params.omega;
%     eps0=params.eps0;
%     mu0=params.mu0;
%     etha0=params.etha0;
%     
%     xlowwidth=pml.xlowwidth;
%     xhighwidth=pml.xhighwidth;
%     ylowwidth=pml.ylowwidth; 
%     yhighwidth=pml.yhighwidth;
%     zlowwidth=pml.zlowwidth;   
%     zhighwidth=pml.zhighwidth;
%     
%     pmltol=1e-8;
%     
%     d=dx+dy+dz;
%     
%     R=pml.R;
%     m=pml.m;
%     sigma.xhigh=-(m+1)*log(R)/(2*etha0*xhighwidth);
%     sigma.xlow=-(m+1)*log(R)/(2*etha0*xlowwidth);
%     sigma.yhigh=-(m+1)*log(R)/(2*etha0*yhighwidth);
%     sigma.ylow=-(m+1)*log(R)/(2*etha0*ylowwidth);
%     sigma.zlow=-(m+1)*log(R)/(2*etha0*zlowwidth);
%     sigma.zhigh=-(m+1)*log(R)/(2*etha0*zhighwidth);
%     
%     pml_xhigh=Lx - xhighwidth;
%     pml_xlow = xlowwidth;
%     pml_yhigh=Ly - yhighwidth;
%     pml_ylow = ylowwidth;
%     pml_zhigh=Lz - zhighwidth;
%     pml_zlow = zlowwidth;
%     
  
%     xhpar=@(x)abs(( x-pml_xhigh)/xhighwidth  );
%     
%     xlpar=@(x)abs(( x-pml_xlow)/xlowwidth  );
%     yhpar=@(y)abs(( y-pml_yhigh)/yhighwidth  );
%     ylpar=@(y)abs(( y-pml_ylow)/ylowwidth  );
   KappaMax=1;
    
        pmltol=params.tol;
    dx=params.dx;    dy=params.dy;    dz=params.dz;
    
    hx=params.hx;    hy=params.hy;    hz=params.hz;

    %     Lx=params.Lx;    Ly=params.Ly;    Lz=params.Lz;
    
    tol=params.tol;
    k0=params.k0;

        
    xlowwidth=pml.xlowwidth;    xhighwidth=pml.xhighwidth;
    ylowwidth=pml.ylowwidth;    yhighwidth=pml.yhighwidth;
    zlowwidth=pml.zlowwidth;    zhighwidth=pml.zhighwidth;
    

    d=dx+dy+dz;

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
    
    xhpar=@(x)abs(( x-pml_xhigh)/xhighwidth  );
    xlpar=@(x)abs(( x-pml_xlow)/xlowwidth  );
    yhpar=@(y)abs(( y-pml_yhigh)/yhighwidth  );
    ylpar=@(y)abs(( y-pml_ylow)/ylowwidth  );
    
    sxf=@(x) ...
              (x      >=  pml_xhigh).*(1 + 1j* xhpar(x).^m * a.xhigh )+...
              (x      <=  pml_xlow ).*(1 + 1j* xlpar(x).^m * a.xlow)+...
              (x      > pml_xlow ).*( x      <pml_xhigh );   
         
%               (xyz(1) >=  pml_xhigh).*(1+(KappaMax-1)*xhpar(xyz(1)) + 1j* xhpar(xyz(1)).^m*sigma.xhigh/(eps0*omega))+...
%                (xyz(1) <=  pml_xlow ).*(1+(KappaMax-1)*xlpar(xyz(1)) + 1j*xlpar(xyz(1)).^m*sigma.xlow /(eps0*omega))+...
%                (xyz(1) > pml_xlow ).*( xyz(1) <pml_xhigh );

        
%                         (xyz(1)>pml_xhigh)*(1+1j*sigmamax/k0*  sin(pi/2*((xyz(1)-pml_xhigh)/xhighwidth))^2)*(1 + (amax-1)*((xyz(1) - pml_xhigh)/ xhighwidth)^p ) +...
%                 (xyz(1)<pml_xlow )*(1-1j*sigmamax/k0*  sin(pi/2*((pml_xlow -xyz(1))/xlowwidth))^2) *(1 + (amax-1)*((xlowwidth - xyz(1))/ xlowwidth )^p ) +...
%                 (xyz(1)>=pml_xlow)*(xyz(1)<=pml_xhigh) ;s
...

     syf=@(xyz) 1;...
%               (xyz(2) >=  pml_yhigh).*(1 + 1j* abs(( xyz(2)-pml_yhigh)/yhighwidth  ).^m*sigma.yhigh/(eps0*omega))+...
%                (xyz(2) <=  pml_ylow ).*(1 + 1j* abs(( xyz(2)-pml_ylow )/ylowwidth   ).^m*sigma.ylow /(eps0*omega))+...
%                (xyz(2) > pml_ylow ).*( xyz(2) <pml_yhigh );
    szf=@(xyz) 1;
%     (xyz(3)>pml_zhigh)*(1-1j*sigmamax*  sin(pi*((xyz(3)-pml_zhigh)/(2*zlowwidth))^2))*(1 + (amax-1)*((xyz(3)  -  pml_zhigh)/zhighwidth)^p ) +...
%                 (xyz(3)<pml_zlow)*(1-1j*sigmamax*  sin(pi*((pml_zlow -xyz(3))/(2*zlowwidth))^2))*(1 + (amax-1)*((zlowwidth-xyz(3))/ zlowwidth)^p ) +...
%                (xyz(3)>=pml_zlow)*(xyz(1)<=pml_zhigh) ;
           
    
%     invsxf=@(xyz) 1/sxf(xyz);
%     invsxf=@(xyz) 1/syf(xyz);
%     invsxf=@(xyz) 1/szf(xyz);
    
    syz_xf=@(xyz) syf(xyz)*szf(xyz)/sxf(xyz);
    sxz_yf=@(xyz) sxf(xyz)*szf(xyz)/syf(xyz);
    sxy_zf=@(xyz) sxf(xyz)*syf(xyz)/szf(xyz);
    
    invsyz_xf=@(xyz) 1/syz_xf(xyz);
    invsxz_yf=@(xyz) 1/sxz_yf(xyz);
    invsxy_zf=@(xyz) 1/sxy_zf(xyz);    
  
% sx=  round(amen_cross(2*ones(d, 1),@(bin) sxf(binarr2xyz(bin-1,...
%      dx,dy,dz,hx,hy,hz)), 1e-14,'nswp',40),tol);  
% sy=  round(amen_cross(2*ones(d, 1),@(bin) syf(binarr2xyz(bin-1,...
%      dx,dy,dz,hx,hy,hz)), 1e-14,'nswp',40),tol); 
% sz=  round(amen_cross(2*ones(d, 1),@(bin) szf(binarr2xyz(bin-1,...
%      dx,dy,dz,hx,hy,hz)), 1e-14,'nswp',40),tol); 
 
syz_x=round(amen_cross(2*ones(d, 1),@(bin) syz_xf(bin2xyz(bin-1,...
     params)), pmltol,'nswp',40),tol);
sxz_y=round(amen_cross(2*ones(d, 1),@(bin) sxz_yf(bin2xyz(bin-1,...
     params)), pmltol,'nswp',40),tol);
sxy_z=round(amen_cross(2*ones(d, 1),@(bin) sxy_zf(bin2xyz(bin-1,...
     params)), pmltol,'nswp',40),tol);
 
invsyz_x=round(amen_cross(2*ones(d, 1),@(bin) invsyz_xf(bin2xyz(bin-1,...
     params)), pmltol,'nswp',40),tol);
invsxz_y=round(amen_cross(2*ones(d, 1),@(bin) invsxz_yf(bin2xyz(bin-1,...
     params)), pmltol,'nswp',40),tol);
invsxy_z=round(amen_cross(2*ones(d, 1),@(bin) invsxy_zf(bin2xyz(bin-1,...
     params)), pmltol,'nswp',40),tol);

end