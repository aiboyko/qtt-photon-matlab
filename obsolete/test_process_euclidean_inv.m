
% 
% A=[12 3 .1; 231 12 1; 3e-2 22 2];
% iA=inv(A);
tol1=1e-4
Xprev=1e-6*A

for i=1:10
i
buf1=round(A*Xprev,tol1)
buf2=round(Xprev*buf1,tol1)
Xnext=round(2*Xprev-buf2,tol1)
Xprev=Xnext
end
%     r(i)=norm(Xprev-iA,1);

% plot(r)

