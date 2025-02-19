%diffraction_grating
% kron(ones(4,1),[0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0])
custom_macromask=kron(ones(4,1),[0 0 1 0]);
phize(ttZ,ones(photonic_crystal_size),jackal_params)
