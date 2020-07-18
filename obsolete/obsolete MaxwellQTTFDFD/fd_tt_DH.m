function y= fd_tt_DH(N,D,isPeriodic,tol,BlochPhase)
%it is implied that N=2;
if D > 1
y=tt_shift(N,D,0)-tt_shift(N,D,1);
    if nargin< 3
    tol=0;
    end
    if isPeriodic
    y=y-BlochPhase*tt_shift(N,D,-(N^D-1));
    end
else
    y=tt_matrix([1 -BlochPhase*isPeriodic;-1 1])
end

y=round(y,tol);
return
end