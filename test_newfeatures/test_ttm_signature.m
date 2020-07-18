AA=tkron(tt_matrix([1; 0]),A)
cA=core2cell(AA)

for i=1:numel(cA)
sz(i,:)=[size(cA{i}) ones(1,4-numel(size(cA{i})))]
end