classdef BagOfPinuses
    properties
        listofpinuses
    end
    methods
        function obj=set.listofpinuses(obj,val)
            keyboard
            'list of pinuses in this bag is being changed'
            obj.listofpinuses=val;
        end 
    end
end