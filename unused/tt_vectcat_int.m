function y = tt_vectcat_ext(varargin)
%By Alexey Boyko, Skolkovo Institute of Science and Technology,
%a.boyko@skoltech.ru
%alexey.boyko@skolkovotech.ru
%Moscow, Russia, 2015

%INPUTS
% arr {1d cell array} of <tt_matrix>  variables
%OUTPUTS
% y <tt_matrix> 

%this function assembles a larger vector in tt-format from a given 1d cell array
%of vectors in tt-format
%so this glues vectors internally 
%this does
% [a1, b1, .., z1,a2,b2,..,z2,..,aN,bN,..zN]
y=tt_reshape(vertcat(varargin{:}),[nargin,size(varargin{1})]); 

return 
end