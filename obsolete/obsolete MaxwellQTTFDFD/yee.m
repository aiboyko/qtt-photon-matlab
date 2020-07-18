function y=yee(f)
    switch f 
        case 'Ex'
            idz=1;
            idy=1;
            idx=2;
        case 'Ey'
            idz=1;
            idy=2;
            idx=1;            
        case 'Ez'
            idz=2;
            idy=1;
            idx=1;              
        case 'Hx'
            idz=2;
            idy=2;
            idx=1;
        case 'Hy'
            idz=2;
            idy=1;
            idx=2;            
        case 'Hz'
            idz=1;
            idy=2;
            idx=2;             
    end
        y= [idz idy idx];
end