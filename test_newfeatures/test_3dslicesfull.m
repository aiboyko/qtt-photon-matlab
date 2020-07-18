for i=1:10:128
    figure;
    imagesc(real(reshape(F3_L_3d(i,:,:),[2^dy,2^dz]).'));
    colorbar;
    axis equal tight
end