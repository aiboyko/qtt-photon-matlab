% dx=4
% dy=dx %important!!
% tol=1e-5;
% fb=rand(2^(dx+dy),1);
% fA=rand(2^(dx+dy));
% 
% b=tt_tensor(reshape(fb,2*ones(1,dx+dy)),tol)
% A=tt_matrix(reshape(fA,2*ones(1,2*(dx+dy))),tol)
% Ab=round(A*b,tol)

list_for_standard_permute=[1:dx;dx+1:dx+dy]; list_for_standard_permute=list_for_standard_permute(:)'
[~,inv_list_for_standard_permute]=sort(list_for_standard_permute) % needed for tt_ overloaded 

% tic
% p_A=round(permute(A,list_for_standard_permute,tol),tol)
% toc
tic
p_G=round(permute(G,list_for_standard_permute,tol),tol)
p_chi=round(permute(chi,list_for_standard_permute,tol),tol)
pp_A=round(tt_eye(2,d)+hx*hy*hz*p_G*diag(p_chi),tol)
toc
% norm(p_A-pp_A)/norm(p_A)

p_b=round(permute(b,list_for_standard_permute,tol),tol)

tic
p_x_qtt=amen_solve2(pp_A,p_b,tol,'resid_damp',1.0,'kickrank',9,'max_full_size',10000,'trunc_norm','residual')
toc
x_qtt2=round(permute(p_x_qtt,inv_list_for_standard_permute,tol),tol)
% pp_ttAb=round(p_ttA*p_ttb,tol)
% ttAv2=permute(pp_ttAb,inv_list_for_standard_permute,tol)
% 
% norm(ttAv2-Ab)/norm(Ab)