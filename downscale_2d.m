function pic_downscaled=downscale_2d(pic)
% This function downscales like this in 2D [ .>][<. ][ .>][<. ]  ->  [   .   ][   .   ]

% pic=[ 1 2 3 4 ;1 2 3 4 ; -1 -2 -3 -4 ;-1 -2 -3 -4]'
N=max(size(pic));
downscale_M=kron(speye(N/2),kron([1 1],kron(speye(N/2),[1 1])))/4;

pic_downscaled=reshape(downscale_M*pic(:),[N/2,N/2]);
end