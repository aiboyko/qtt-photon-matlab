wz=wavesol.x;
wz=tt_subtensor(wavesol.x,d+1,CONST.AXIS.Z);
ttv2=tt_subtensor(wz,1:(dz+dy),[ones(1,dz) 1+de2bi(ytest,dy)]);
figure;plot(abs(full(ttv2) ))