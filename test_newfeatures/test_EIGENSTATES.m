opts.isreal=0;
opts.tol=1e-8;
opts.disp=0;
opts.maxit=50;
opts.isreal=0;
%opts.p=200;
%opts.p
NPLx=4;
NPLy=4;
[V,D,eigflag] = eig(Axh,2^d,NPLx*NPLy,'li',opts);
eigenval=diag(D);
figure;
for i=1:NPLx
    for j=1:NPLy
        
        counter=j+NPLy*(i-1);
        subplot(NPLx,NPLy,counter)
        imagesc(real(reshape(V(:,counter),[2^dx 2^dy])));axis equal tight;colormap(jet(2^12));
        %colorbar
        title(['eigv=',num2str(eigenval(counter))])
    end
end
