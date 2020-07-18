function y=downscale_2d_qtt(v,params)
% This function downscales like this in 2D [ .>][<. ][ .>][<. ]  ->  [   .   ][   .   ]

% pic=[ 1 2 3 4 ;1 2 3 4 ; -1 -2 -3 -4 ;-1 -2 -3 -4]'
dx=v.d*params.dx/params.d;%what is thst supposed to mean? î_Î
dy=v.d*params.dy/params.d; 

% dx=params.dx;
% dy=params.dy;

tol=params.tol;

ttm11=tt_matrix([1 1],tol);

downscale_M=mtkron(   ttm11,tt_eye(2,dy-1),ttm11,tt_eye(2,dx-1) )/4; %not sure about dx and dy order here =( Needs further check

y=purify(round(downscale_M*v,1e-14)); %purify strips the ttv from 1x1 modes, that appeared after rectangular matvec

end