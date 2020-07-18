function y=Dforw(N,h)
aux_m_diag=diag(ones([1, N]),0)
aux_m_diag_for=diag(ones([1, N]),1);
aux_m_diag_for=aux_m_diag_for(1:N,1:N);
y=(aux_m_diag_for-aux_m_diag)/h;
end