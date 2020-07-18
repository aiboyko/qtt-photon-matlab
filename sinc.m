function y=sinc(x)

y=zeros(size(x))
    y(x == 0)=1;
    xgood=x(x ~= 0);
    y(x ~= 0)=sin(xgood)./(xgood);
end