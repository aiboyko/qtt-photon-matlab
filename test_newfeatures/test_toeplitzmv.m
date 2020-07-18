clear all

N=1000;
c1=[rand(N,1)];
c2=[c1(1);rand(N-1,1)];
T=toeplitz(c1,c2);
tc=[c1;flipud( c2(2:numel(c2)))];
x=rand(N,1);

xext=[x;zeros(N-1,1)];
mvC=@(c,x)ifft(fft(c) .*fft(x) );
keyboard 

tic
fullMV=T*x;
toc
tic
tempMV=mvC(tc,xext);
tempMV= tempMV(1:N);
toc

norm(fullMV-tempMV)

% BTTB entities
% main entity - 2d table of numbers (2N-1)(2M-1)
% 1) full via custom function
% 2) fast via fft


A= rand(5,1);
A=kron(ones(1,5),A);
for i=1:5
    A(:,i)=i*A(:,i);
    A(1,i)=11*i;
end

BTTB=blocktoeplitz(A);
x=rand(9,1);

mv=BTTB*x;

z=zeros(5);
z(1:3,1:3)=reshape(x,[3 3]);
mvfast=ifft2(fft2(A).*fft2(z));
mvfast=reshape(mvfast(1:3,1:3),[9,1]);

norm(mvfast-mv)


