%% --Amen-preconded GMRES block
% % 
% disp('holy15')
% holy_precond15=@(x)full(amen_solve2(A,tt_tensor(reshape(x,2*ones(1,d)),tol),tol,'nswp',15,'kickrank',4,'verb',1 ,'resid_damp',1.0));
% tic
% [x1_holy15,flag,relres_holy15,iter,resvec_holy15]=gmres(Axh,ffb,100,tol,100,holy_precond15);
% toc
% figure;plot(log10(resvec_holy15/norm(ffb)),'-o');
% xlabel('iter');
% ylabel('log10 of relres_holy15');

% % % % % holy=@(x)holy_precond(A,x,d,tol,36);
% % % % % tic
% % % % % [x1_holy10,flag,relres_holy5,iter,resvec_holy5]=gmres(Axh,...
% % % % %     ffb,5,tol,5,holy);
% % % % % toc
% % % % % figure;plot(log10(resvec_holy5/norm(ffb)),'-o');
% % % % % xlabel('iter');
% % % % % ylabel('log10 of relres_holy');

% % % % % figure('name','GMRES-holy FFT solution');imagesc(abs(reshape(  x1_holy10  ,[Nx,Ny])));colorbar;axis equal tight;
% % % % % colormap(cool(1024)); 
% % % % % caxis([0,3])
% % % % % hold on
% % % % % % contour(abs(reshape(full(chi/k0^2/eps),[Nx,Ny])),'LineWidth',2,'Color','w')
% % % % % % contour(abs(reshape(full(source),[Nx,Ny]))*1e-10,'Color','r')
% % % % % hold off

% % 
% % disp('holy3')
% % holy_precond3=@(x)full(amen_solve2(A,tt_tensor(reshape(x,2*ones(1,d)),tol),tol,'nswp',3,'kickrank',5,'verb',0 ));
% % tic
% % [x1_holy2,flag,relres_holy3,iter,resvec_holy3]=gmres(Axh,ffb,10,tol,10,holy_precond3)
% % toc
% % figure;plot(log10(resvec_holy3/norm(ffb)),'-o');
% % xlabel('iter');
% % ylabel('log10 of relres_holy3');