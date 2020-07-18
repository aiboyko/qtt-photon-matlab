%this bydlocode produces a Fermi-like circle on pixelized geometry such
%that integral of the pixel circle is equal to pi*R^2 up to circtol
clear all;
close all;
filename = 'circleOO.gif'
circtol=1e-12;


MhFermi=@(r,R,kT)1./ ( 1+exp(( r-R )/kT )); %double(E-EF <1.5*kT)*
    XLEFT=-2;
    XRIGHT=2;
    YLEFT=-2;
    YRIGHT=2;

    Lx=XRIGHT-XLEFT; 
    Ly=YRIGHT-YLEFT;

dx0=2;
dxmax=10;
imax=dxmax-dx0+1;
for i=1:imax
    i
    R=Lx*0.25;
    kT=R*7e-2/(2^(i-1));
    dx=dx0-1+i;
    dy=dx;
    d=dx+dy;

    Nx=2^dx; 
    Ny=2^dy;

    hx=Lx/(Nx);
    hy=Ly/(Ny);

    % |[ . ][ . ][ . ][ . ]|
    % |<--------L--------->|
    %
    Nx=2^dx; 
    Ny=2^dy;

    Xcenter=(XRIGHT+XLEFT)/2;
    Ycenter=(YRIGHT+YLEFT)/2;

    [XX,YY]=meshgrid(XLEFT+ hx/2+ (0:(2^dx-1))*hx,YLEFT+ hy/2+(0:(2^dy-1))*hy);
    %___Geometry___
    Rs=sqrt((XX-Xcenter).^2+(YY-Ycenter).^2);

    % if i==1
    %     R=Lx*0.25;
    %     kT=R*1e-3;
    % else
    %     R=Lx*0.25*sqrt(oversize(i-1));
    %     kT=R*1e-3*sqrt(oversize(i-1));
    % end

    %if a circle is 32x32pix, the coef should be 1e-2
    relerrS(i)=1;
    Sth=pi*R^2;
    count=0;
    while (relerrS(i)>circtol)&&count<500
        count=count+1;
        ffmask=MhFermi(Rs,R,kT);
   
        S=sum(sum(ffmask))*hx*hy;
        oversize(i)=abs(S)/abs(Sth);
        relerrS(i)=abs(S-Sth)/abs(Sth);
        R=R*sqrt(1/oversize(i));
        kT=kT*sqrt(1/oversize(i));
    end
    relerrS(i)=abs(S-Sth)/abs(Sth);
    figure();imagesc(XX(1,:),YY(:,1),ffmask);axis equal tight;caxis([0,1]); colorbar; colormap(gray(2^12)); 
%     drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
      if i == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end

end

figure;plot(dx0:dxmax,log10(relerrS),'-o');ylabel('log10(relerrS)');xlabel('dx');