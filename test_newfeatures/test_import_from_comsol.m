fname='G:\comsoldata.txt'
% tx=textread(fname,'%s','commentstyle','matlab');

% dat=cellfun(@(x)str2num(x),tx);
% dat=reshape(dat,[N^2,3]);
% X=(reshape(dat(:,1),[N,N]));
% Y=(reshape(dat(:,2),[N,N]));
% F=(reshape(dat(:,3),[N,N]));

dat=dlmread(fname);
N=sqrt(numel(dat)/3);
dat=reshape(dat,[3,numel(dat)/3]);

dat=reshape(dat,[N^2,3]);
X=(reshape(dat(:,1),[N,N]));
Y=(reshape(dat(:,2),[N,N]));
F=(reshape(dat(:,3),[N,N]));

figure;(imagesc(abs(F)'));axis equal tight;colormap(jet(2^10));caxis([0,2])