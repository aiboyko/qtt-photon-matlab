function ttGXY=abinitioG_entity(params,varargin)
%This function creates matrix G from Helmholtz Volume Integral equation
%from first principles

        hx=params.grid.hx;
        hy=params.grid.hy;
        dx=params.grid.dx;
        dy=params.grid.dy;
        dz=params.grid.dz;
        mej=params.quadratures.maxerrorjumps;
        d=dx+dy+dz;
        k0=params.rhs.k0;
        tol=params.tol;
        softeningcoef=params.quadratures.softeningcoef;

        if dz==0
        ttdeltaX=hx*(tt_x(2,dx+1)-2^dx);%volume-centric lattice is assumed [ . ][ . ]
        ttdeltaY=hy*(tt_x(2,dy+1)-2^dy);

        %general approach with composition of functions
%       gh_curv2xy=@(coord) [coord(:,1).*cos(coord(:,3)),coord(:,2).*cos(coord(:,4)), coord(:,1).*sin(coord(:,3)) , coord(:,2).*sin(coord(:,4))];

        %direct computation with cosine theorem
        %gh_total=@(coord) 1./sqrt( coord(:,1).^2+coord(:,2).^2 - 2*coord(:,1).*coord(:,2).*cos(coord(:,3)-coord(:,4))  +   tol*1e1)

        %generating standard 4D mesh
        ttdeltaXmesh=mtkron(ttdeltaX,tt_ones(2,dy+1)); %correct!
        ttdeltaYmesh=mtkron(tt_ones(2,dx+1),ttdeltaY);
        
        if nargin==1
%             gh=@(xy)hx*hy*gh_xy_Nyst_selfheal(xy,params);
             gh=@(xy)hx*hy*gh_xy_Nyst_shifted(xy,params);
        else
            keyboard
            %USELESS! Cross can't adapt for sufficiently singular patch, it
            %just reeplaces it with large but different number even on high
            %precisions 
            LCN_patch=varargin{1};
            gh=@(xy) hx*hy*gh_xy_LCN(xy,LCN_patch,params);
        end

        ttGXY=round(amen_cross({ttdeltaXmesh,ttdeltaYmesh},@(x)gh(x),tol*1e-2,...
        'trunc_method','svd','kickrank',100, 'nswp',100,'max_err_jumps',mej,'zrank',20),tol);
        %figure('name','ttmatrix with custom reshape and LCN function handle');imagesc(abs(full(G)));colormap(jet(2^12));colorbar;axis equal tight
    else
        error('3D problems are not yet supported')
    end
    
end