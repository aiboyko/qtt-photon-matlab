function y=fake_hankel(s,L,k0)
%this is an analytic formula for fourier image of truncated Hankel
% fake_hankel creates stuff with centered (!!) zero
% therefore you need to use fftshift before usage of fft

% which wavevectors should I sample?


y= (1+...
    1j*pi/2*L * s.* besselj(1,L*s) * besselh(0,1,L*k0)-...
    1j*pi/2*L * k0* besselj(0,L*s) * besselh(1,1,L*k0))./(s.^2 - k0^2);
end

