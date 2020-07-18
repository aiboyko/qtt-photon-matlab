%jackel compression. custom TT
I = imread('jackel2.png');
tol=1e-5;
Z=zeros(256);
Z(1:size(I,1),1:size(I,2))=double((255-I(:,:,1)))/255;
ttZ=tt_tensor(reshape(Z,2*ones(1,log2(numel(Z)))),tol)
figure;imagesc(reshape(full(Z),[256 256]));axis equal tight;caxis([0 16]);colormap(jet(1024))

piss=zeros(4);
piss(3,1)=1;
tt_piss=tt_tensor(reshape(piss, 2*ones(1,4) ),1e-5);


ttZZ=tkron(tt_piss,ttZ)

figure;imagesc(reshape(full(ttZZ),[2^10,2^10]));caxis([0 16]);colormap(jet(1024));axis equal tight

figure;imagesc(reshape(full(tt_subtensor(ttZZ,[9 10 19 20],[ 2 2 1 1])),[2^8,2^8]));caxis([0 16]);colormap(jet(1024));axis equal tight

% figure;imagesc(reshape(ZZ_abomination(:,1,2,1,1),[2^8 2^8]));caxis([0 16]);colormap(jet(1024));axis equal tight
% 
ZZ_abomination=reshape(ZZ,[size(Z,1)*size(Z,2) 2*ones(1,Zdx)  2*ones(1,Zdx)]);
ttZZ2=tt_tensor(ZZ_abomination,tol)
figure;imagesc(reshape(full(tt_subtensor(ttZZ2,[3 4 5],[ 1 1 1])),[2^10,2^7]));caxis([0 16]);colormap(jet(1024));axis equal tight

% figure;imagesc(ZZ);axis equal tight;caxis([0 16])
% % colormap(flipud(gray))
% colormap(jet(1024))
