% In seems natural to me that timeseries should be recorded
% as a horzcat of tt_tensor. 
% One must hope that spatial similarity of the field state in different time
% slices will contribute to keeping the total tensor rank low

%UPD nope, mem grows linearly, and time feels like quadratically

%UPD cell mem grows also linearly, but time is also linear ->profit.

tol=1e-14;
ts=round(eps_Ex,tol);
td=12;
mmemz=zeros([65 1]);
memz(1)=mem(eps_Ex);
% for i=1:2^td-1
%     ts=horzcat(ts,i*sin(2*pi*i/256)/100*eps_Ex);
%     ts=round(ts,tol);
% %     memz(i+1)=mem(ts);
% end
% ts
% mem(ts)
% ts2=reshape(ts,2*ones([1 8+td]));
% ts2=round(ts2,1e-14)
% mem(ts2)
c={};

for i=1:100000;
    c{i}=eps_Ex*i;
end