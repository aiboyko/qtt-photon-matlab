function pic_upscaled=upscale_2d(pic)
% This function upscales like this in 2D     [  = . =  ][ = . = ]  -> [ .>][<. ][ .>][<. ]

% pic=[ 1 2 3 4 ;1 2 3 4 ; -1 -2 -3 -4 ;-1 -2 -3 -4]'
% N=max(size(pic));
% downscale_M=kron(speye(N/2),kron([1 1],kron(speye(N/2),[1 1])))/4;

pic_upscaled=kron(pic,[1 1;1 1]); %use fro matrix norm in case you dont convert this to vector

end