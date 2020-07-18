function [ s ] = gaussian_source( t,t0,sigma )
%GAUSSIAN Summary of this function goes here
%   Detailed explanation goes here

s=exp(-(t-t0).^2/(2*sigma^2));




end

