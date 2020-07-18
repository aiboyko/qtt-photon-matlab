%jackel compression. custom TT
I = imread('jackel64_2.png');
tol=1e-6;
Z=zeros(64);
jackal_params.dx=log2(size(Z,2));
jackal_params.dy=log2(size(Z,1));
jackal_params.tol=1e-4;

photonic_crystal_size=2^9;

Z(1:size(I,1),1:size(I,2))=double((255-I(:,:,1)))/255;
Z=Z*numel(Z)/sum(sum(Z));

[XX,YY]=meshgrid(-photonic_crystal_size/2:photonic_crystal_size/2-1,  -photonic_crystal_size/2:photonic_crystal_size/2-1);
sigma=photonic_crystal_size/4;
gausss=exp(-sqrt(XX.^2 + YY.^2)/sigma);
figure;imagesc(gausss)
 axis equal tight;
 caxis([0 1]);
 colormap(magma(1024))

ttZ=tt_tensor(reshape(Z,2*ones(1,log2(numel(Z)))),tol);


figure;
imagesc(reshape(full(Z),2.^[jackal_params.dx jackal_params.dy]));
 colormap(magma(1024))
 axis equal tight;
 colorbar
 
ttj2=phize(ttZ,gausss,jackal_params)
%ttj3=tt_subtensor(ttj2,[1:6 ,17:17+5],ones(1,12))

% ttj3=tt_subtensor(ttj2,[1:6 ,ttj2.d/2:ttj2.d/2+5],ones(1,12))

ttj3=tt_sumovermodes(ttj2,[1:6 ,ttj2.d/2+1:ttj2.d/2+6])/(2^12)

fullJackals=full(ttj3);
 figure;imagesc(reshape(fullJackals,sqrt(numel(fullJackals))*[ 1 1]));
 axis equal tight;caxis([0 1]);colormap(magma(1024))
%  caxis([0 1])
% 
% piss=zeros(4);
% piss(3,1)=1;
% tt_piss=tt_tensor(reshape(piss, 2*ones(1,4) ),1e-5);
% 
% ttZZ=tkron(tt_piss,ttZ)
% 
% figure;imagesc(reshape(full(ttZZ),[2^10,2^10]));caxis([0 16]);colormap(jet(1024));axis equal tight
% 
% figure;imagesc(reshape(full(tt_subtensor(ttZZ,[9 10 19 20],[ 2 2 1 1])),[2^8,2^8]));caxis([0 16]);colormap(jet(1024));axis equal tight
% 
% % figure;imagesc(reshape(ZZ_abomination(:,1,2,1,1),[2^8 2^8]));caxis([0 16]);colormap(jet(1024));axis equal tight
% % 
% ZZ=full(ttZZ);
% ZZ_abomination=reshape(ZZ,[size(Z,1)*size(Z,2) 2*ones(1,Zdx)  2*ones(1,Zdx)]);
% ttZZ2=tt_tensor(ZZ_abomination,tol)
% figure;imagesc(reshape(full(tt_subtensor(ttZZ2,[3 4 5],[ 1 1 1])),[2^10,2^7]));caxis([0 16]);colormap(jet(1024));axis equal tight

% figure;imagesc(ZZ);axis equal tight;caxis([0 16])
% % colormap(flipud(gray))
% colormap(jet(1024))
