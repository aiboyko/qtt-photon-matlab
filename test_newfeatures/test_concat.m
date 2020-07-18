N=16
% ttv00=tt_tensor((1:N).^(-1),1e-16,[2 2 2 2 ]);
ttv1=tt_tensor(1:8,1e-16,[2 2 2 ]);
ttv2=tt_tensor((1:16).^2,1e-16, [2 2 2 2]);
ttv=tt_reshape(horzcat(ttv1,ttv2),[2 2 2 2 2 ]);
%vertcat can be reshaped only if you add

ttv_TP2=tt_reshape(vertcat(ttv0,ttv1,ttv2),[3 2 2 2 2 ]); %this does
% [a1, b1, .., z1,a2,b2,..,z2,..,aN,bN,..zN]


ttv_TP=tt_reshape(horzcat(ttv0,ttv1,ttv2),[2 2 2  2 3]); %so this glue vectors externally 
%aka [a1,..,aN,b1,..,bN,..,z1,..zN]

ttv=tt_vectcat_ext(ttv0,ttv1,ttv2)