% field_data=wavesol_exp.x{numel(wavesol_exp.x)};
field_data=wavesol_exp.x{4}
dxmax=8;
dymax=8;
fun=@abs;
%% extracting and plotting Electric Field - z
[ttfield,newdx,newdy,newdz ]=qtt_Esubfield(field_data,params,CONST.AXIS.Z,dxmax,dymax,0);
[X,Y]=meshgrid(0:2^newdy-1,0:2^newdx-1);
fff=reshape(full(ttfield),[2^newdy 2^newdx]);
fff=permute(fff, [2 1] );
figure;

imagesc(0:2^newdy-1,0:2^newdx-1,fun(fff));colorbar
axis equal tight
hold on
%% extracting and plotting material

if eps~=1
    ttmat=qtt_subfield(mask_of_coords,params,dxmax,dymax,0);
    mmm=reshape(full(ttmat),[2^newdy 2^newdx]);
    mmm=permute(mmm, [2 1] );
    contour(X,Y,real(mmm),[1-.1 1+.1],'-w','LineWidth',2);
    axis equal tight
end
%%
if ~doUPML
ttsx=qtt_subfield(invsx,params,dxmax,dymax,0);
else
ttsx=qtt_subfield(invsxy_z,params,dxmax,dymax,0);
end
sss=reshape(full(ttsx),[2^newdy 2^newdx]);
sss=permute(sss, [2 1] );
sss=imag(1./sss);
contour(X,Y,real(sss),[sqrt(tol) sqrt(tol)],'-.g','LineWidth',2)
axis equal tight
%%
ttq=qtt_subfield(box_Q,params,dxmax,dymax,0);
qqq=reshape(full(ttq),[2^newdy 2^newdx]);
qqq=permute(qqq, [2 1] );
contour(X,Y,real(qqq),[sqrt(tol) sqrt(tol)],'-.r','LineWidth',2)
axis equal tight

hold off
caxis([0,2])