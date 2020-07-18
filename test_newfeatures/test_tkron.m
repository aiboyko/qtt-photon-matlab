%this file is to test by-hand assembly of Laplace operator via
%Kronecker product.

d=3;
ttl=tt_qlaplace_dd([d]);
ttl2=tt_qlaplace_dd([d d]);
ttl2_artificial= tkron(ttl,tt_eye(2,d)) + tkron (tt_eye(2,d), ttl);
% norm(full(ttl2_artificial)-full(ttl2))

ttl3=tt_qlaplace_dd([d d d]);

ttl3_artificial=mtkron(ttl,tt_eye(2,d),tt_eye(2,d)) +...
mtkron(tt_eye(2,d), ttl,tt_eye(2,d))+...
mtkron(tt_eye(2,d), tt_eye(2,d),ttl);

% norm(full(ttl3_artificial)-full(ttl3))

%example of piece-wise tensor assembly 
%tkron 
art=tkron(ttl,tt_matrix([ 1 0;0 0]))
