ttr=diag(tt_rand(2,3,1))
for i=1:8
    for j=1:8
        ttr_elwise(i,j)=tt_submatrix(ttr,[1:3],indexify([i j],ttm_coresizes(ttr,'only_sizes')) );
    end
end
norm(full(ttr)-ttr_elwise)

