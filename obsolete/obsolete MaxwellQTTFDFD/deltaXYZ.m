function [ deltaXYZ ] = deltaXYZ( field )
%DELTAXYZ this is a shifts vector for Yee staggered grid
       switch(field)
       case 'Hx' 
           deltaXYZ=[.5 0 0];
       case 'Hy' 
           deltaXYZ=[0 .5 0];
       case 'Hz' 
           deltaXYZ=[0 0 .5];
       case 'Ex'
           deltaXYZ=[0 .5 .5];
       case 'Ey' 
           deltaXYZ=[.5 0 .5];
       case 'Ez' 
           deltaXYZ=[.5 .5 0];  
       end

end

