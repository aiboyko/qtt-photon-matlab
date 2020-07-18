function y=PMLcoefs(params)
%test prototype (unfinished)
%lossless and AMEnless PML constructor

low=10
high=10
m=2
pml_low=round(round(1- tt_heaviside(2,dx,low+1),tol).* round((-tt_x(2,dx)+low)/low,tol).^m,tol)
pml_high=round(round(tt_heaviside(2,dx,2^dx-high+1),tol) .* round( ( (tt_x(2,dx)-2^dx+high+1)/high).^m,tol),tol)
pmlovski=round(pml_low+pml_high,tol)
figure;plot(full(tt_elem_reverse(1+pmlovski)))
end