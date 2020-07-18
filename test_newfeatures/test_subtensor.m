
d=10;
y=tt_sin_cos(d,2*pi/(2^d),1)

figure(1);plot(full(y))
G=core2cell(y);
G

list_of_idx = [1:6];
list_of_vals = ones(1,4); %should besame size as list_of_idx

Nidx=numel(list_of_idx)
%this thing freezes 7th harmonic to 1st value (not the 2nd), and calculates
%subtensor, which is smaller

for k=Nidx:-1:1
    idx = list_of_idx(k)
    val = list_of_vals(k)
    BackConvolve = (idx>1); 
    G{idx}=G{idx}(:,val,:)
    
    if numel(G{idx})~=1
        if BackConvolve
            %here is the scenario fo convolve with previous index   
            G{idx}=permute(G{idx},[1 3 2])
            sz=size(G{idx-1});
            G{idx-1}=reshape(G{idx-1},[sz(1)*sz(2),sz(3)]);
            G{idx-1}=G{idx-1}*G{idx};
            G=G(setdiff(1:numel(G), idx) ); %we keep all core numbers except for the frozen one (idx)
            G{idx-1}=reshape(G{idx-1},[sz(1),sz(2),size(G{idx-1},2)]);
        else
            %here is the scenario fo convolve with forward index           
            G{idx}=permute(G{idx},[1 3 2]);
            sz=size(G{idx+1});
            if numel(sz)==2
                sz(3)=1;
            end
            G{idx+1}=reshape(G{idx+1},[sz(1),sz(2)*sz(3)]);
            G{idx+1}=G{idx}*G{idx+1};
            G=G(setdiff(1:numel(G), idx) );
            %we keep all core numbers except for the frozen one (idx)
            G{idx}=reshape(G{idx},[size(G{idx},1),sz(2),sz(3)]);
        end
    end
        
end
if numel(G) ==1
    y1=G{:};
else
    y1=cell2core(tt_tensor,G);
end

figure;plot(full(y1),'-o')
