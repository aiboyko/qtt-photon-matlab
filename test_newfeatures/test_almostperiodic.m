%% creating macromask
 mdx=3;
 mdy=3;
% macromask=ones(2^mdx,2^mdy);
% macromask(2^(mdx-1),:)=0


macromask = [1 1 1 1 1 1 1 1;...
             0 0 0 0 0 0 0 1;...
             1 1 1 1 1 1 0 1;...
             1 1 1 1 1 1 0 1;...
             0 0 0 0 1 1 0 1;...
             1 1 1 0 1 1 0 1;...
             1 1 1 0 1 1 0 1;...
             1 1 1 0 1 1 0 1];
    

%% applying macromsk
tmacromask=tt_tensor(reshape(macromask,2*ones(1,mdx+mdy)),tol);
material_mask = tt_tensor(reshape(ffmask,2*ones(1,d)),tol);
tmacromask=round(tkron(material_mask,tmacromask),tol)
newmask=round(permute(tmacromask,[1:dx (d+1):(d+mdx) (dx+1:d) (d+mdx+1):(d+mdx+mdy)],tol),tol)
figure; imagesc(reshape(full(newmask),[2^(dx+mdx),2^(dy+mdy)]));axis equal tight