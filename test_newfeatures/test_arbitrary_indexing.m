function y=indexify(xval,list_of_coresizes)
    for i=1:size(list_of_coresizes,2)
        size_x=list_of_coresizes(1,i);
        idx_x(i)=mod(xval,size_x);
        xval=xval-idx_x(i);
        xval=xval/size_x;
    end
    y=idx_x;
end