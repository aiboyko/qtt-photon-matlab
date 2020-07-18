classdef classTTBinaryMask < classTTField
    %making a class that has tt-structures as its fields
    %this 
        
    methods
        function obj=classTTBinaryMask(dx,newttv)
            obj=obj@classTTField(dx,newttv);
        end
        function obj=and(a,b)
            obj=a;
            obj.ttv=and(a.ttv,b.ttv);
        end
        function obj=or(a,b)
            obj=a;
            obj.ttv=or(a.ttv,b.ttv);
        end
        function obj=not(a)
            obj=a;
            obj.ttv=not(a.ttv);
        end
        function obj=xor(a,b)
            obj=a;
            obj.ttv=xor(a.ttv,b.ttv);
        end
        function visualizeField(obj)
            visualizeField@classTTField(obj);
            colormap(gray(2^10));
            caxis([0,1]);
            axis equal tight;
        end
    end
    
end

