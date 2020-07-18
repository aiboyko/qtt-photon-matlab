function [ r ] = test_r( xyz )
%TEST_R Summary of this function goes here
%   Detailed explanation goes here
r=sqrt(xyz(:,1).^2+xyz(:,2).^2+xyz(:,3).^2);

end

