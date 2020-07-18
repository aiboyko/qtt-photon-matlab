function res= helm3d_kernel(s,k0,L)

    div_zero_idx_k0 = abs(s-k0)<1e-5;
    div_zero_idx_0 = abs(s)<1e-5;
    s(div_zero_idx_k0) = k0 + 1e-5;
    s(div_zero_idx_0) = 1e-5;
    
    res = (-1 + exp(1.0j*L*k0).*(cos(L*s)-1.0j*k0./s.*sin(L*s)))./(k0^2-s.^2);
    
    lim_re = -1.0/(2*k0) * (-L*cos(L*k0)*sin(L*k0) - (sin(L*k0))^2/k0 + L*(cos(L*k0))^2);
    lim_im = -1.0/(2*k0) * (-L + sin(L*k0)*cos(L*k0)/k0);
    
    res(div_zero_idx_k0) = lim_re + 1.0j*lim_im;
    res(div_zero_idx_0) = (-1.0 + exp(1.0j*L*k0)*(1-1.0j*k0*L))/k0^2;
    
end