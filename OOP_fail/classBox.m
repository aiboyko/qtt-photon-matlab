classdef classBox < handle % so if the box is a handle it's convenient to make shared boxes, 
    %but it's unclear how to detect box change to recalc the Part tt and CS tt
    properties
        boundaries
    end
    properties(SetAccess=protected)
        ISUP2DATE
        myTTBinaryMask
    end
    methods
        
        function obj=classBox(val)
            obj.boundaries=val;
        end
        function y=boundaries.set(obj,val)
            if ~isequal(val,obj.boundaries)
            	y.boundaries=val;
                y.myTTBinaryMask=[];
            end
        end
        function y=myTTBinaryMask.get(obj)
            if ~isempty(obj.myTTBinaryMask)
                %SK, please calc my myTTBinaryMask
                %return myTTBinaryMask
                %1) callback?
                %2) params? and calculate here
            else
                y=obj.myTTBinaryMask;
            end
        end        
    end
end