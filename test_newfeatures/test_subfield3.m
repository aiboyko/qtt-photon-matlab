

wz=tt_subtensor(wavesol.x,d+1,CONST.AXIS.Z);

ttv_crossection_refl=tt_subtensor(wz,dz+dy+1:dz+dy+dx ,[1+de2bi(xtest,dx)]);
crossection_refl=full(ttv_crossection_refl);

ttv_crossection_trans=tt_subtensor(wz,dz+dy+1:dz+dy+dx ,[1+de2bi(xtest_trans,dx)]);
crossection_trans=full(ttv_crossection_trans);

fcross_refl=reshape(crossection_refl,[2^dz 2^dy]);
fcross_trans=reshape(crossection_trans,[2^dz 2^dy]);

R=sum(abs(crossection_refl ).^2)/numel(crossection_refl)
T=sum(abs(crossection_trans).^2)/numel(crossection_trans)
E=T+R
log_err=log10(abs(1-E))

% figure
% imagesc(abs(fcross_refl));
% colorbar
% figure
% imagesc(abs(fcross_trans));
% colorbar