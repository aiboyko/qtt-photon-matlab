function [ v_reshaped,X,Y,Z ] = F1_3DtoF1_L( v3d )
%QTT2FIELD Summary of this function goes here

v_reshaped=permute(v3d,[3 2 1]);
v_reshaped=reshape(v_reshaped,[numel(v3d) 1]);

end
