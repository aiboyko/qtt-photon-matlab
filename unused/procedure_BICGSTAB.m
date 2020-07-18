%% -BICGSTAB solution block
%     disp('bicgstab pure')
%     tic
%     [x,flag,relres,iter,resvec2]=bicgstab(Axh,ffb,tol,500);
%     toc
%     figure;plot(log10(resvec2/norm(ffb)),'-o');
%     xlabel('iter');
%     ylabel('log10 of relres2');
%     data2=real(reshape(  x  ,[Nx,Ny,Nz]));
%     data2=data1(:,:,2^(dz/2));
%     figure('name','Bicgstab FFT solution');imagesc(data2);colorbar;axis equal tight;
%     colormap(jet(1024)); 
%     caxis([-1,1])
% %     hold on
% %     contour(abs(reshape(full(chi/k0^2/eps),[Nx,Ny,Nz])),'LineWidth',2,'Color','w')
% %     % contour(abs(reshape(full(source),[Nx,Ny]))*1e-10,'Color','r')
% %     hold off
% 
%     keyboard