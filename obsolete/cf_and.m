function [ y ] = cf_and( xyz, C )
%ALL_IDIOTS_UNITE Summary of this function goes here
%   Detailed explanation goes here
y=1;
for i=1:size(C,2)

    f=C{i};
    y=y&&f(xyz);
end

end

