function y=ttm_zeros(N,D)
%By Alexey Boyko, Skolkovo Institute of Science and Technology,
%a.boyko@skoltech.ru
%alexey.boyko@skolkovotech.ru
%Moscow, Russia, 2015

%returns a zero matrix N^D by N^D in tt_matrix format

y=diag(tt_zeros(N,D));
return
end