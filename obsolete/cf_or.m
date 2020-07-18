function [ y ] = cf_or( xyz, C )
%ALL_IDIOTS_UNITE Summary of this function goes here
%   Detailed explanation goes here
y=0;
for i=1:size(C,2)

    f=C{i};
    y=y||f(xyz);
end

end

