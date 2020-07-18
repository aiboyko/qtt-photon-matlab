classdef classCoordinateSystem
%basic class for coordinate system
%carries information of some structured grid
%what is the range in units of length
%units of length is NANOMETRES
    
    properties
        ds %[dx dy]
        bordervalues %[xleft, xright; yleft, yright ]
        boxes %list of boxes
        parts %list of details
    end
    
    
    properties (Dependent)
        N %Total numbers of grid elements towards each dimension
        L %Total lengths of the domain in each dimension
    end
% %     
    properties (Access=protected)

        boxes_ttmaskvault %list of tt_vectors with recpective binary masks. SHould always correspond to the list of boxes
        parts_ttmaskvault
    end
    
    methods
        function obj=classCoordinateSystem(bordervalues,ds)
            obj.bordervalues=bordervalues;            %[xleft, xright; yleft, yright ]
            obj.ds=ds; %[dx dy]
        end
        function value=get.L(obj)
            value=(obj.bordervalues(:,2)-obj.bordervalues(:,1)).';
        end
        function value=get.N(obj)
            value=2.^obj.ds;
        end 
   
    end
    
end

