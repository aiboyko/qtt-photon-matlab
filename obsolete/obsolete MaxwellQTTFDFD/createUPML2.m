function [syz_x,sxz_y,sxy_z,invsxz_y,invsyz_x,invsxy_z]=createUPML2(params,pml)
%unfinished
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
     KappaMax=1
%     
%     xhpar=@(x)abs(( x-pml_xhigh)/xhighwidth  );
%     xlpar=@(x)abs(( x-pml_xlow)/xlowwidth  );
%     yhpar=@(y)abs(( y-pml_yhigh)/yhighwidth  );
%     ylpar=@(y)abs(( y-pml_ylow)/ylowwidth  );
%     
%     ktxhigh=sigma.xhigh/(eps0*omega)
%     ktxlow=sigma.xlow /(eps0*omega)
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
%     sigma.xhigh=-(m+1)*log(R)/(2*etha0*xhighwidth)
%     sigma.xlow=-(m+1)*log(R)/(2*etha0*xlowwidth)
%     sigma.yhigh=-(m+1)*log(R)/(2*etha0*yhighwidth);
%     sigma.ylow=-(m+1)*log(R)/(2*etha0*ylowwidth);
%     sigma.zlow=-(m+1)*log(R)/(2*etha0*zlowwidth);
%     sigma.zhigh=-(m+1)*log(R)/(2*etha0*zhighwidth);
%     
%     switch params.units
%     case 'SI'
%         ktxhigh=sigma.xhigh/(eps0*omega)
%         ktxlow=sigma.xlow/(eps0*omega) 
%     case 'Gauss'
%         ktxhigh=sigma.xhigh/(k0)
%         ktxlow=sigma.xlow/(k0)     
%     end

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
    
    sxf=@(x)...
               (x      >=  pml_xhigh).*(1 + 1j* xhpar(x).^m * a.xhigh )+...
               (x      <=  pml_xlow ).*(1 + 1j* xlpar(x).^m * a.xlow)+...
               (x      > pml_xlow ).*( x      <pml_xhigh );    
    
%     sxf=@(x) ...
%               (x >=  pml_xhigh).*(1+(KappaMax-1)*xhpar(x) + 1j* xhpar(x).^m*ktxhigh)+...
%                (x <=  pml_xlow ).*(1+(KappaMax-1)*xlpar(x) + 1j*xlpar(x).^m*ktxlow)+...
%                (x > pml_xlow ).*( x  <pml_xhigh );
           
%     syf=@(y) 1;
%     szf=@(z) 1;
    
    invsxf=@(x) 1/(sxf(x)+1e-9 );
    
sx=amen_cross(2*ones(dx, 1),@(bin) sxf( bi2de(bin-1)*hx ), pmltol,'nswp',40,'trunc_method','cross')
inv_sx=amen_cross(2*ones(dx, 1),@(bin) invsxf( bi2de(bin-1)*hx ), pmltol,'nswp',40)

%     sy=  round(amen_cross(2*ones(d, 1),@(bin) syf(binarr2xyz(bin-1,...
%      dx,dy,dz,hx,hy,hz)), 1e-14,'nswp',40),tol); 
%     sz=  round(amen_cross(2*ones(d, 1),@(bin) szf(binarr2xyz(bin-1,...
%      dx,dy,dz,hx,hy,hz)), 1e-14,'nswp',40),tol); 

    
    
% %we will use hard-coded X-only PML    
%     syz_xf=@(xyz) 1/(sxf(xyz)+pmltol);
%     sxz_yf=@(xyz) sxf(xyz);
%     sxy_zf=@(xyz) sxf(xyz);
%     
%     %sxf(xyz)=mtkron(tt_ones(2,dz),tt_ones(2,dy),sxf
%     
%     invsyz_xf=@(xyz) sxf(xyz)/(syf(xyz)*szf(xyz)+pmltol);
%     invsxz_yf=@(xyz) syf(xyz)/(sxf(xyz)*szf(xyz)+pmltol);
%     invsxy_zf=@(xyz) szf(xyz)/(sxf(xyz)*syf(xyz)+pmltol);   


%     syz_xf=@(xyz) syf(xyz)*szf(xyz)/(sxf(xyz)+pmltol);
%     sxz_yf=@(xyz) sxf(xyz)*szf(xyz)/(syf(xyz)+pmltol);
%     sxy_zf=@(xyz) sxf(xyz)*syf(xyz)/(szf(xyz)+pmltol);
%     
%     invsyz_xf=@(xyz) sxf(xyz)/(syf(xyz)*szf(xyz)+pmltol);
%     invsxz_yf=@(xyz) syf(xyz)/(sxf(xyz)*szf(xyz)+pmltol);
%     invsxy_zf=@(xyz) szf(xyz)/(sxf(xyz)*syf(xyz)+pmltol);    
  
 
% syz_x=round(amen_cross(2*ones(d, 1),@(bin) syz_xf(bin2xyz(bin-1,...
%      params)), pmltol,'nswp',40,'trunc_method','cross'),tol);
% sxz_y=round(amen_cross(2*ones(d, 1),@(bin) sxz_yf(bin2xyz(bin-1,...
%      params)), pmltol,'nswp',40,'trunc_method','cross'),tol);
% sxy_z=round(amen_cross(2*ones(d, 1),@(bin) sxy_zf(bin2xyz(bin-1,...
%      params)), pmltol,'nswp',40,'trunc_method','cross'),tol);
% keyboard
syz_x=purify(mtkron(tt_ones(2,dz),tt_ones(2,dy),inv_sx));
sxz_y=purify(mtkron(tt_ones(2,dz),tt_ones(2,dy),sx));
sxy_z=purify(mtkron(tt_ones(2,dz),tt_ones(2,dy),sx));


% invsyz_x=round(amen_cross(2*ones(d, 1),@(bin) invsyz_xf(bin2xyz(bin-1,...
%      params)), pmltol,'nswp',40,'trunc_method','cross'),tol);
% invsxz_y=round(amen_cross(2*ones(d, 1),@(bin) invsxz_yf(bin2xyz(bin-1,...
%      params)), pmltol,'nswp',40,'trunc_method','cross'),tol);
% invsxy_z=round(amen_cross(2*ones(d, 1),@(bin) invsxy_zf(bin2xyz(bin-1,...
%      params)), pmltol,'nswp',40,'trunc_method','cross'),tol);

invsyz_x=purify(mtkron(tt_ones(2,dz),tt_ones(2,dy),sx));
invsxz_y=purify(mtkron(tt_ones(2,dz),tt_ones(2,dy),inv_sx));
invsxy_z=purify(mtkron(tt_ones(2,dz),tt_ones(2,dy),inv_sx));
end