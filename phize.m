function newmask=phize(material_mask,macromask,params)
% material_mask is tt_tensor
% macromask is a MATLAB binary array 8x8
%%
%this 
 mdx=log2(size(macromask,1));
 mdy=log2(size(macromask,2));
 tol=params.tol;
 dx=params.dx;
 dy=params.dy;
 d=dx+dy;
 
 
% macromask=ones(2^mdx,2^mdy);
% macromask(2^(mdx-1),:)=0
if ~isequal(class(macromask),'tt_tensor')
    tmacromask=tt_tensor(reshape(macromask,2*ones(1,mdx+mdy)),tol);
else
    tmacromask-macromask
end
% material_mask = tt_tensor(reshape(ffmask,2*ones(1,d)),tol);
newmask=round(tkron(material_mask,tmacromask),tol);
newmask=round(permute(newmask,[1:dx (d+1):(d+mdx) (dx+1:d) (d+mdx+1):(d+mdx+mdy)],tol),tol)
% figure; imagesc(reshape(full(newmask),[2^(dx+mdx),2^(dy+mdy)]));axis equal tight
end