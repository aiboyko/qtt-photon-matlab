function y=tt_lininterp2(ttv,tol,params)
%this function takes a tt_vector and creates a 2 times longer

%1 2 3 4 | 5 6 7 8 9 | 10 11 12 13 14 15| 16 17
%(1) 2 3 4 5 | (6) 7 8 9 10 11 | (12) 13 14 15 16 17 | 18 19
    dz=params.dz;
    dy=params.dy;
    dx=params.dx;


    dd=ttv.d;
 
%     ttv=tt_sin_cos(dd,2*pi/(2^dd),1);
    G=core2cell(ttv);
    G2=cell([dd+3 1]);
    
    G2(1)={1};
    G2(2:dz+1)=G(1:dz);
    
    G2(dz+2)={1};
    G2(dz+3:dz+dy+2)=G(dz+1:dz+dy);
    
    G2(dz+dy+3)={1};
    G2(dz+dy+4:dz+dy+dx+3)=G(dz+dy+1:dz+dy+dx); 
    
    G2(dz+dy+dx+4:dd+3)=G(dz+dy+dx+1:dd); 
    
    ttv2=cell2core(tt_tensor,G2);

    coreup=tt_matrix([1;0],tol);
    coredown=tt_matrix([0;1],tol);
   

    one_eyed=tt_matrix([1 0;0 0],tol);
    list={};
    
    ds=[dx dy dz]
    for j=1:3
        lowerdiag=tt_shift(2,ds(j),1);
        for i=1:ds(j)
            list{i}=one_eyed;
        end

        big_one_eyed=round(mtkron(list),tol);

        interpmat1d{1}=round(0.5*tkron(coreup,big_one_eyed)+0.5*(tkron(coreup,tt_eye(2,ds(j)))+tkron(coreup,lowerdiag)) +...
            tkron(coredown,tt_eye(2,ds(j))),tol);
    end
%     full(interpmat);
keyboard
     y=round(interpmat*ttv2,tol);
%     figure;plot(full(ttv2),'-x');
end