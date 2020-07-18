%test_rectangle_fft

%rectangle func is: 
%0 outside of |x|<=1/2
%1 inside
%analytically,fourier of it is sinc(pi x) = sin (pi x ) / (pi x)

xleft=-.5; %use -2lambda, +2 lambda for comparison with Mie
L=1;
xright=xleft+L;


