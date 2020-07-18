function y = tt_vectcat_ext(varargin)
%By Alexey Boyko, Skolkovo Institute of Science and Technology,
%a.boyko@skoltech.ru
%alexey.boyko@skolkovotech.ru
%Moscow, Russia, 2015

%INPUTS
% arr {1d cell array} of <tt_matrix>  variables
%OUTPUTS
% y <tt_matrix> 

%! vectors should be equal sizes
%this function assembles a larger vector in tt-format from a given 1d cell array
%of vectors in tt-format
%so this glues vectors externally 
%aka [a1,..,aN,b1,..,bN,..,z1,..zN]

y=tt_reshape(horzcat(varargin{:}),[size(varargin{1}), nargin]); 

return 
end