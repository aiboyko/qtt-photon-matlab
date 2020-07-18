classdef classHelloHandle < handle
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        veshalka
    end
    
    properties (Access=private)
        private_veshalka
    end
    
    properties (Access=protected)
        protected_veshalka
    end
    
    methods
        function obj=classHelloHandle(val)
            obj.veshalka=val;
            obj.protected_veshalka=val^2;
            obj.private_veshalka=val^3;
        end
        function changePriV(obj,val)
            obj.private_veshalka=val
        end
        
        function changeProV(obj,val)
            obj.protected_veshalka=val
        end
        
        function changeV(obj,val)
            obj.veshalka=val
        end
    end
    
end

