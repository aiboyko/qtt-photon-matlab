classdef BagOfHandlePinuses
    properties
        listofhandlepinuses
    end
    methods
        function obj=set.listofhandlepinuses(obj,val)
            
            'list of pinuses in this BagOfHandlePinuses is being changed'
            obj.listofhandlepinuses=val;
        end 
    end
end