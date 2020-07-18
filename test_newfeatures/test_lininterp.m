
dd=3;
tol=1e-6;
ttv=tt_sin_cos(dd,2*pi/(2^dd),0);
X1=(0:2^dd-1)/(2^dd-1)*(2^(dd+1)-1);
figure;plot(X1,full(ttv),'-bo');
hold on
     
    dd=ttv.d;
    G=core2cell(ttv);
    G(2:dd+1)=G;
    G(1)={1};
    
    
    ttv2=cell2core(tt_tensor,G);

    coreup=tt_matrix([1;0],tol);
    coredown=tt_matrix([0;1],tol);
    lowerdiag=tt_shift(2,dd,1);

    one_eyed=tt_matrix([1 0;0 0],tol);
    list={};
    for i=1:dd
        list{i}=one_eyed;
    end
    big_one_eyed=round(mtkron(list),tol);
%0.5*tkron(coreup,big_one_eyed)+
    interpmat=round(0.5*(tkron(coreup,tt_eye(2,dd))+tkron(coreup,lowerdiag)) +...
        tkron(coredown,tt_eye(2,dd)),tol);
%     full(interpmat);

     y=round(interpmat*ttv2,tol);
     X2=(0:2^(dd+1)-1);
     plot(X2,full(y),'-rx');
     hold off