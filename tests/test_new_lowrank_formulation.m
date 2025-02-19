alpha=1e-5;
newA=round(tt_blockwise({tt_eye(2,params.grid.d),G/alpha;alpha*diag(chi),-tt_eye(2,params.grid.d)}),params.tol)
keyboard
newb=tkron(b,tt_tensor([1 0],1e-12,[2]));
tic
new_x_qtt=amen_solve2(newA,newb,params.tol,'resid_damp',1.2,'kickrank',12,'max_full_size',15000,'trunc_norm','residual','local_prec','rjacobi');
toc

pic_qtt=abs(reshape(full(new_x_qtt),[2^params.grid.dx,2^params.grid.dy]));
pic_qtt=permute(pic_qtt,[2 1]);




figure('name','new Amen solution');
ffmask=reshape(full(material_mask),...
    [params.grid.Ny,params.grid.Nx]);
if exist('XX')
    imagesc(XX(1,:),YY(:,1).',pic_qtt);
    hold on
    contour(XX,YY,ffmask.',1,'LineWidth',2,'Color','k')
    contour(XX,YY,ffmask.',1,'LineWidth',1,'Color','w')
    hold off
else
    imagesc(pic_qtt);   
    hold on
    contour(ffmask.',1,'LineWidth',2,'Color','k')
    contour(ffmask.',1,'LineWidth',1,'Color','w')
    hold off
end
colorbar;
colormap(magma(2000))
caxis(2*[-0 1])
axis equal tight
ylabel('x');
xlabel('y');
title('new qtt') 




f=@(x)1./(x(:)); 
divine=amen_cross({chi},f,params.tol*1e-1,'trunc_method','svd','kickrank',100, 'nswp',100,'max_err_jumps',15,'zrank',20)