d1 = 2;
d2 = 2;
N1 = 2 ^ d1;
N2 = 2 ^ d2;

 a = reshape([1:4*N1*N2], [2*N1, 2*N2])
%a=[1:4*N1*N2];
aq = tt_tensor(reshape(a, 2*ones(1, d1 + d2 + 2)));
af = tt_qtoepl(aq, [d1, d2]);
round(full(af))