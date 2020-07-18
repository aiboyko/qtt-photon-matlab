function piece = tt_vectsplit_ext(ttv,k)
%By Alexey Boyko, Skolkovo Institute of Science and Technology,
%a.boyko@skoltech.ru
%alexey.boyko@skolkovotech.ru
%Moscow, Russia, 2015

%INPUTS
% ttv <tt_tensor> 2^d * M
%OUTPUTS
% y <tt_tensor> 

%this function extracts k-th part of the original vector

%what it was: [v1_1,..,v1_(2^d),v2_1,..,v2_(2^d),..,vM_1,..vM_(2^d)] <-tt_tensor
%what it became: [vk_1,...,vk_N]  <-tt_tensor

lastCoreIdx=max(size(size(ttv)));
piece=ttv;
piece{lastCoreIdx}=piece{lastCoreIdx}(:,k);

return 
end