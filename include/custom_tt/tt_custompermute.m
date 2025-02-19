function ttv_custom=tt_custompermute(ttv,tol)

%dxdydxdy
%ex1
% 1 2 3 4 5 6 | 7 8 || ... | ... -> 1 7 2 8 3 4 5 6 ||
% 1 2 3 | 4 5 6 7 8 || ... | ... -> 1 4 2 5 3 6 7 8 ||
    

     d=ttv.d/2;
     
%     dm=min(dx,dy)
%     
%     list1=1:dx;
%     list2=dx+1:dx+dy;
%     list_for_custom_permute1=[list1(1:dm);list2(1:dm)];
%     list_for_custom_permute1=list_for_custom_permute1(:)';
%     list_for_custom_permute1=[list_for_custom_permute1,list1(dm+1:dx),list2(dm+1:dy)];
%     
%     list3=dx+dy+1:2*dx+dy;
%     list4=2*dx+dy+1:2*dx+2*dy;    
%     list_for_custom_permute2=[list3(1:dm);list4(1:dm)];
%     list_for_custom_permute2=list_for_custom_permute2(:)';
%     list_for_custom_permute2=[list_for_custom_permute2,list3(dm+1:dx),list4(dm+1:dy)];    
    
    % list_for_custom_permute=[list_for_custom_permute1,list_for_custom_permute2]
    %     list_for_custom_permute1=[list1(1:dm);list2(1:dm)];

    list_for_standard_permute=[1:d;d+1:2*d];
    list_for_standard_permute=list_for_standard_permute(:)';
%     list_for_custom_permute=[1     6     2     7     3     8     4     9     5    10]
    ttv_custom=permute(ttv,list_for_standard_permute,tol);
end