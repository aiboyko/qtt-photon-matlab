classdef HandlePinus < handle
    properties
        length
    end
    methods
        
        function obj=set.length(obj,val)
            'length of this HandlePinus is being changed'
            obj.length=val;
        end 
        function y=get.length(obj)
            'length of this HandlePinus is being read'
            y=obj.length;
        end 
        function [ sumpinus ] = plus( pinus1,pinus2 )
            sumpinus.length=pinus1.length+pinus2.length;
        end
    end
end