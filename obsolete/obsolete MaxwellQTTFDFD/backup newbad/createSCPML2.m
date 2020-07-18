function [invsx,invsy,invsz]=createSCPML2(params,pml)
%Attention! coordinates translation here is done manually, via *hx, change
%for staggered grid

    pmltol=1e-10;
    dx=params.dx;    dy=params.dy;    dz=params.dz;
    
    hx=params.hx;    hy=params.hy;    hz=params.hz;
    
%     Lx=params.Lx;    Ly=params.Ly;    Lz=params.Lz;
    
    tol=params.tol;
    k0=params.k0;
    omega=params.omega;
    eps0=params.eps0;
    mu0=params.mu0;
    etha0=params.etha0;
        
    xlowwidth=pml.xlowwidth;    xhighwidth=pml.xhighwidth;
    ylowwidth=pml.ylowwidth;    yhighwidth=pml.yhighwidth;
    zlowwidth=pml.zlowwidth;    zhighwidth=pml.zhighwidth;
    
%     pml_xhigh=Lx - xhighwidth;
%     pml_xlow = xlowwidth;
%     pml_yhigh=Ly - yhighwidth;
%     pml_ylow = ylowwidth;
%     pml_zhigh=Lz - zhighwidth;
%     pml_zlow = zlowwidth;
    
    xhpar=@(x)abs(( x-pml_xhigh)/xhighwidth  );
    xlpar=@(x)abs(( x-pml_xlow)/xlowwidth  );
    yhpar=@(y)abs(( y-pml_yhigh)/yhighwidth  );
    ylpar=@(y)abs(( y-pml_ylow)/ylowwidth  );

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

    
    sxf2=@(x) ...
              (x      >=  pml_xhigh).*(1 + 1j* xhpar(x).^m * ktxhigh )+...
              (x      <=  pml_xlow ).*(1 + 1j* xlpar(x).^m * ktxlow)+...
              (x      > pml_xlow ).*( x      <pml_xhigh );           

    syf2=@(y) 1;
    szf2=@(z) 1;

    invsxf2=@(x) 1/sxf2(x);
    invsyf2=@(y) 1/syf2(y);
    invszf2=@(z) 1/szf2(z);
    

invsx=  round(mtkron(tt_ones(2,dz),tt_ones(2,dy),amen_cross(2*ones(dx, 1),@(bin) invsxf2( bi2de(bin-1)*hx ), pmltol,'nswp',40)),tol);
invsy= round(mtkron(tt_ones(2,dz),amen_cross(2*ones(dy, 1),@(bin) invsyf2(bi2de(bin-1)*hy), pmltol,'nswp',40),tt_ones(2,dx)),tol);
invsz= round(mtkron(      amen_cross(2*ones(dz, 1),@(bin) invszf2(bi2de(bin-1)*hz), pmltol,'nswp',40)   ,tt_ones(2,dy),tt_ones(2,dx)  ),tol);
% 
% keyboard

       invsx_full=F1_LtoF1_3D(full(invsx), params);
%        invsx_full2=F1_LtoF1_3D(full(invsx2), params);
%        
%        figure();imagesc(real( invsx_full(:,:,1)));colorbar; 
%        figure();imagesc(imag( invsx_full(:,:,1)));colorbar;
%        
%        figure();imagesc(real( invsx_full2(:,:,1)));colorbar; 
%        figure();imagesc(imag( invsx_full2(:,:,1)));colorbar;
%        
       figure;plot(real( invsx_full(:,:,1)))
       figure;plot(imag( invsx_full(:,:,1)))
%        figure;plot(real( invsx_full2(:,:,1)))
%        figure;plot(imag( invsx_full2(:,:,1)))
       
keyboard
      
end