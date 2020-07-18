% tts6=round(tt_x(2,6).*tt_x(2,6),1e-8)

d=6;
tts6=tt_sin_cos(d,2*pi/(2^d),1)
tts6=tts6.*tts6
% figure(1);plot(full(tts6),'-ro')
G=core2cell(tts6)

% G1=G{1}
% G2=G{2}
G(2:7)=G
G1=[1 1]
% G1=G1(:,:,2)
% G2=G2(2,:,:)
G(1)={G1}
% G(2)={G2}


figure(2);plot(full(cell2core(tt_tensor,G)),'-o')
