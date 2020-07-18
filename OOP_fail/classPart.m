classdef classPart
    % in fact this class is a prototype of an SK(CS) class, since boxes are not
    % supposed to be stored inside parts.
    % SK(CS) does part-level mask assembly, realclassPart itself only deletes its part_ttfield_mask when it's obsolete 
    % SK(CS) does box-level mask assembly,  realclassBox  itself only deletes its box_ttfield_mask when it's obsolete 
    properties
        boxes %classBox is supposed to be here, but for now it's made for just classTTField. However, if those are handles
        composition;
        eps;
        %logical string like 'boxmask{1}&~boxmask{2}' implied
    end
    
    properties (SetAccess=protected)
        overall_TTBinaryMask
    end
    
    methods
        function y=get.overall_TTBinaryMask(obj)
            if isempty(obj.boxes)
                'no boxes to make the overall_mask'
                y=[];

            elseif isempty(obj.composition)
                'no composition to make overall_mask with'
                y=[];
            elseif isempty(obj.overall_TTBinaryMask)
%                 keyboard
                eval(strcat('composer=@(boxmask)' ,obj.composition,';'));
                obj.overall_TTBinaryMask=composer(obj.boxes);
                y=classTTBinaryMask(2,    round(obj.overall_TTBinaryMask.ttv,5e-14)); %kinda system constant, ~5e-14 for TT-Toolbox
            else
                y=obj.overall_TTBinaryMask;
            end

        end
        function obj=set.boxes(obj,val)
            obj.boxes=val;
            
            if ~isempty(obj.overall_TTBinaryMask)
                obj.overall_TTBinaryMask=[];          
            end
        end        
        
    end
    
end