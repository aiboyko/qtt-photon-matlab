classdef classTTField<handle
    %making a class that has tt-structures as its fields
    % BETTER TO RELOAD ENTIRE BOOL OPERATORS ON tt_tensor
    properties
        ttv
        dx
    end
    properties (Dependent)
        dy
        d
    end    
    
    methods
        %constructor
        function obj= classTTField(dx,newttv)
            obj.dx=dx;            
            obj.ttv=newttv;
        end
        %dependent numbers
        function buf=get.dy(obj)
            buf=obj.d-obj.dx;
        end   
        function buf=get.d(obj)
            buf=obj.ttv.d;
        end
        %methods
        
        function visualizeField(obj,varargin)
            if nargin==1
                figure();
            else
                figure(varargin{1})
            end
            fttv=reshape(full(obj.ttv),2.^[obj.dx,obj.dy]);
            imagesc(real(fttv));
            colormap(jet(2^10));
            axis equal tight; 
            colorbar;
%             caxis([-100,100]);
        end
        %overloaded operators
%         function c=and(a,b)
%             c.dx=a.dx;
%             c.ttv=and(a.ttv,b.ttv);
%         end
%         function c=or(a,b)
%             c.dx=a.dx;
%             c.ttv=or(a.ttv,b.ttv);
%         end
%         function c=not(a)
%             c.dx=a.dx;
%             c.ttv=not(a.ttv);
%         end
%         function c=xor(a,b)
%             c.dx=a.dx;
%             c.ttv=xor(a.ttv,b.ttv);
%         end

    end
    
end

