
n=15
nremove=10;

%it appears that it's absolutely mandatory to measure axis from 0, not from
%1, otherwise even is everything is perfect, plots will be shifted

sinful=tt_sin_cos(n,2*pi*2^-(n-1),1)
figure;
hold on
plot([0:2^(n)-1],full(sinful),'-gx')


sincore=core2cell(sinful)
newsincore=sincore(nremove+1:n)

core1=newsincore{1}
core11=core1(1,:,:)

newsincore(1)={core11}

% newsincore2=newsincore;
% core12=core1(1,:,:)
% newsincore2(1)={core12}
newsinful = cell2core(tt_zeros(2,n),newsincore)



plot((2^nremove)*[0:2^(n-nremove)-1],full(newsinful)','-bo')
hold off

%% 3D decompression
%vector of EH format
%d=15, 7+7+1
%extraction of different field components
enum_field=1;
enum_component=3;
c=core2cell(sol.x)
bufc_field=c{d+2}(:,enum_field,:)
c(d+2)=[]
bufc_component=c{d+1}(:,enum_component,:)
c(d+1)=[]
s=size(bufc_component)
bufc_component=reshape(bufc_component,[prod(s(1:numel(s)-1)),s(numel(s))])*bufc_field
bufc_last=c{d}
s=size(bufc_last);

bufc_last=reshape(bufc_last,[prod(s(1:numel(s)-1)),s(numel(s))])*bufc_component
bufc_last=reshape(bufc_last,s(1:numel(s)-1))
c{d}=bufc_last
ttc=cell2core(tt_zeros(2,d),c)


