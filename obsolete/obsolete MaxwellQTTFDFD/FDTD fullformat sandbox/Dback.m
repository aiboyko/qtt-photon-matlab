function y=Dback(N,h)
aux_m_diag=diag(ones([1, N]),0);
aux_m_diag_back=diag(ones([1, N]),-1);
aux_m_diag_back=aux_m_diag_back(1:N,1:N);
y=(-aux_m_diag_back+aux_m_diag)/h;

end