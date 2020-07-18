function y=tt_getelement(ttv,i)
%% Is redundant
%% may be done with ttv(bin_i)
%example: y=tt_x(2,4)
%bin_i = [1 2 1 1]
    c=core2cell(ttv);
    bufc=c;
    bufprod=1;
    n=numel(c);
    bin_i=de2bi(i,n)+1;

    for k=1:n
        bufcore=bufc{k};
        bufcore=bufcore(:,bin_i(k),:);
        bufc{k}=permute(bufcore,[1 3 2]); %to compress with the neighbour it has to be converted to a matrix
        bufprod=bufprod*bufc{k};
    end
y=bufprod;    

end