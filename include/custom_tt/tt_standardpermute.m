function ttv_custom=tt_standardpermute(ttv,tol)
    d=ttv.d/2;
    list_for_standard_permute=[1:d;d+1:2*d];
    list_for_standard_permute=list_for_standard_permute(:)'

%  Here is what works for converting STANDARD ORDERING 
% note:   
%   list_for_standard_permute=[ 1     2     3     4   [1     2     3     4] + 4]
    ttv_custom=permute(ttv,list_for_standard_permute,tol);
end