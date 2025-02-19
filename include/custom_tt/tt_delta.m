function y=tt_delta(n,L,ind)
y=round(tt_heaviside(n,L,ind)-tt_heaviside(n,L,ind+1),1e-14)
end