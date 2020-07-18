function tt_vec2mat
d=2 %in real life this would be dx+dy+dz, since we are forming matrix d*d
%important!!
tol=1e-5;
x=tt_x(2,d+d);
xONE=reshape(full(x),[2^d,2^d])

list_for_standard_permute=[1:d;d+1:d+d];
list_for_standard_permute=list_for_standard_permute(:)';
x_readyformat=reshape(permute(x,list_for_standard_permute,tol),4*ones(1,d));
ttm_x=tt_matrix(x_readyformat,2*ones(1,d),2*ones(1,d));

xTWO=full(ttm_x)

norm(xONE-xTWO)


% %Permute(ttv)->fullmat
% list_for_probe_permute=[dx+dy:-1:1];
% xALT=permute(x,list_for_probe_permute,tol)%additional permutation
% xALTONE=reshape(full(xALT),[2^dx,2^dy])
% list_for_standard_permute=[1:dx;dx+1:dx+dy]; 
% %ttv->|ttm->permute(ttm)->fullmat
% list_for_standard_permute=list_for_standard_permute(:)';
% ttm_xALT2=permute(ttm_x,list_for_probe_permute,tol)
% xALTTWO=full(ttm_xALT2)
end