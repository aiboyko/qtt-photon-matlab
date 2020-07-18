function G=abinitioG(params,varargin)
%This function creates matrix G from Helmholtz Volume Integral equation
%from first principles

        hx=params.grid.hx;
        hy=params.grid.hy;
        dx=params.grid.dx;
        dy=params.grid.dy;
        dz=params.grid.dz;
        d=dx+dy+dz;
        k0=params.rhs.k0;
        tol=params.tol;

        if dz==0
        ttX=hx/2+hx*tt_x(2,dx);%volume-centric lattice is assumed
        ttY=hy/2+hy*tt_x(2,dy);

        %general approach with composition of functions

%          gh_curv2xy=@(coord) [coord(:,1).*cos(coord(:,3)),coord(:,2).*cos(coord(:,4)), coord(:,1).*sin(coord(:,3)) , coord(:,2).*sin(coord(:,4))]

        %direct computation with cosine theorem
        %gh_total=@(coord) 1./sqrt( coord(:,1).^2+coord(:,2).^2 - 2*coord(:,1).*coord(:,2).*cos(coord(:,3)-coord(:,4))  +   tol*1e1)

        %generating standard 4D mesh
        ttX1=mtkron(ttX,tt_ones(2,dy),tt_ones(2,dx),tt_ones(2,dy)); %correct!
        ttY1=mtkron(tt_ones(2,dx),ttY,tt_ones(2,dx),tt_ones(2,dy));
        ttX2=mtkron(tt_ones(2,dx),tt_ones(2,dy),ttX,tt_ones(2,dy));
        ttY2=mtkron(tt_ones(2,dx),tt_ones(2,dx),tt_ones(2,dy),ttY);

        %doing magical permute
        pttX1=tt_standardpermute(ttX1,tol);
        pttX2=tt_standardpermute(ttX2,tol);
        pttY1=tt_standardpermute(ttY1,tol);
        pttY2=tt_standardpermute(ttY2,tol);

        % curvilinear case
        % pttv=round(amen_cross({pttX1,pttX2,pttY1,pttY2},@(x)gh(gh_curv2xy(x)),tol,...
        % 'trunc_method','svd','kickrank',100, 'nswp',50,'max_err_jumps',3,'zrank',10),tol);
        
%       old
%         if nargin==1
%             gh=@(xxyy) hx*hy*gh_xxyy_Nyst_selfheal(xxyy,params);
%         else
%             LCN_patch=varargin{1};
%             gh=@(xxyy) gh_xxyy_LCN(xxyy,LCN_patch,params);             
%         end

        %new. selfhealing-only
        gh=@(xxyy) hx*hy*gh_xxyy_Nyst_selfheal(xxyy,params);

        pttv=round(amen_cross({pttX1,pttX2,pttY1,pttY2},@(x)gh(x),tol*1e-1,...
        'trunc_method','svd','kickrank',100, 'nswp',50,'max_err_jumps',4,'zrank',10),tol);
        modes = pttv.n(1:d);
        G=round(tt_vec2mat(pttv),tol);
        %figure('name','ttmatrix with custom reshape and LCN function handle');imagesc(abs(full(G)));colormap(jet(2^12));colorbar;axis equal tight
 
    else
        error('3D problems are not yet supported')
    end
    
end