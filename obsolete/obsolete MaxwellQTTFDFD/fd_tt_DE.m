function  y = fd_tt_DE(N,D,isPeriodic,tol,BlochPhase)
%it is implied that N=2;
if D > 1
    y=-tt_shift(N,D,0)+tt_shift(N,D,-1);
    if isPeriodic
    y=y+BlochPhase*tt_shift(N,D,(N^D-1));
    end
else
    y=tt_matrix([-1 1;BlochPhase*isPeriodic -1])
end


y=round(y,tol);
return
end
