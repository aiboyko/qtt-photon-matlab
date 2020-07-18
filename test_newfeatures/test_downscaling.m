data=[ 1 2 3 4 ;1 2 3 4 ; -1 -2 -3 -4 ;-1 -2 -3 -4]'
N=size(data,1);

downscale_M=kron(eye(N/2),kron([1 1],kron(eye(N/2),[1 1])))/4;
reshape(downscale_M*data(:),[N/2,N/2])