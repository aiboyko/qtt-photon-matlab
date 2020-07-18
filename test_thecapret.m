% test_thecapret
macromask= [1 1 1 1; 1 0 0 1; 1 0 0 1; 1 1 1 1]
tol=1e-10;
 mdx=log2(size(macromask,1));
 mdy=log2(size(macromask,2));
tmacromask=tt_tensor(reshape(macromask,2*ones(1,4)),tol);

m=1;
tt_mask=tmacromask;
d=tt_mask.d;
dx=tt_mask.d/2;
dy=tt_mask.d/2;
for k=2:5
    
% m=kron(m0,m);

tt_mask=round(tkron(tt_mask,tmacromask),tol);

tt_mask=round(permute(tt_mask,[1:dx (d+1):(d+mdx) (dx+1:d) (d+mdx+1):(d+mdx+mdy)],tol),tol)
d=tt_mask.d;
dx=tt_mask.d/2;
dy=tt_mask.d/2;
% tt_tensor(reshape(m, 2*ones(1,log2(numel(m)))))


end

% dx1 dx2 dy1 dy2 dx1 dx2 dy1 dy2
% 1   2   3    4   5   6   7   8
% 1 2 5 6 3 4 7 8
% 
% dx1 dx2 dx3 dx4 dy1 dy2 dy3 dy4 dx5 dx6 dy5 dy6
% 1    2   3   4   5   6   7   8   9  10  11  12
% 
% dx1 dx2 dx3 dx4 dy1 dy2 dy3 dy4 dx5 dx6 dy5 dy6
% 1    2   3   4 9 10  5   6   7   8   9  10  11  12

