function BTTB=blocktoeplitz(A)
%create a block Toeplitz matrix with Toeplitz blocks
% A is (2N-1)*(2M-1) array of data, which is exactly what is needed to
% assemble BTTB matrix with MxM block-Toeplitz  matrix with NxN-sized
% Toeplitz blocks

N=(size(A,1)+1)/2;
M=(size(A,2)+1)/2;

BTTB=zeros(M*N);

for b=1:2*M-1
        a=A(:,b);
        c1=a(1:N);
        c2=[a(1); flip(a([N+1:2*N-1]))];
        toeplitz(c1,c2);
        if b<=M
            BTTB=BTTB+kron( diag(ones(M-abs(b-1),1),-b+1) , toeplitz(c1,c2));
        else
            BTTB=BTTB+kron( diag(ones(M-abs(2*M-b),1),2*M-b) , toeplitz(c1,c2));    
        end
end

