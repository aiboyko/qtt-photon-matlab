function y=holy_precond(A,x,d,tol,rmax)
coef=1e-2;
ttx=round(tt_tensor(reshape(x,2*ones(1,d)),0),0,rmax);


tty=amen_solve2(A,ttx,tol*coef,...
    'nswp',1,'kickrank',0,'verb',1,...
    'x0',round(ttx,0,rmax)+tol*tt_ones(2,d),0) ;

y=full(tty);
end