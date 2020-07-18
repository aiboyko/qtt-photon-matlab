fG = full(G);
v = zeros(32,32);
v(1,32) = 1;
v = reshape(v, 1024,1)

figure(1)
mv = fG * v;
mv2 = reshape(mv, 32,32);
imagesc(real(mv2))
colorbar()

figure(2)
mv = transpose(fG) * v;
mv2 = reshape(mv, 32,32);
imagesc(real(mv2))
colorbar()