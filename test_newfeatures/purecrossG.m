function y=purecrossG()
    %Purecross implemented successfully!
    clear all
    dx=5;
    dy=dx;
    Nx=2^dx;
    Ny=2^dy;
    tol=1e-7;

    Xcenter=Nx/2;
    Ycenter=Ny/2;

    ttX=(tt_x(2,dx)-Xcenter)/(Nx-1);
    ttY=(tt_x(2,dy)-Ycenter)/(Ny-1);

    gh=@(x) 1./sqrt(( (x(:,1)-x(:,3)).^2+(x(:,2)-x(:,4)).^2+tol*1e1))

    ttX1=mtkron(tt_ones(2,dy),tt_ones(2,dx),tt_ones(2,dy),ttX); %correct!
    ttY1=mtkron(tt_ones(2,dy),tt_ones(2,dx),ttY,tt_ones(2,dx));
    ttX2=mtkron(tt_ones(2,dy),ttX,tt_ones(2,dy),tt_ones(2,dx));
    ttY2=mtkron(ttY,tt_ones(2,dx),tt_ones(2,dy),tt_ones(2,dx));

    p_ttX1=tt_custompermute(ttX1,dx,dy,tol);
    p_ttX2=tt_custompermute(ttX2,dx,dy,tol);
    p_ttY1=tt_custompermute(ttY1,dx,dy,tol);
    p_ttY2=tt_custompermute(ttY2,dx,dy,tol);
    pttv=amen_cross({p_ttX1,p_ttX2,p_ttY1,p_ttY2},gh,tol,...
    'trunc_method','svd','kickrank',100, 'nswp',70,'max_err_jumps',30,'zrank',10)
    pttv=round(pttv,tol)


    d=pttv.d/2; 
    modes = pttv.n(1:d);
    ttm_custom=round(tt_matrix(reshape(pttv,modes.^2),modes,modes),1e-5)

    t=tt_vec2mat(honest_ttv_custom,tol)
    figure('name','ttmatrix with custom reshape');imagesc(full(ttm_custom));colormap(jet(2^12));colorbar;axis equal tight
    caxis([0,6])
end