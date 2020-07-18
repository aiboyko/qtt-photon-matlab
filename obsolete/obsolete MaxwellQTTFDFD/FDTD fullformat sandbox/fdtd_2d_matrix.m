% DEx
% DEy=Dforw(Ny,hy);
% DEz=Dforw(Nz,hz);
% 
% Db(1,N)=-letmebePeriodic;
% Df(N,1)=letmebePeriodic;
% 
% C_Ex=DEy*Ez - DEz*Ey;
% C_Ey=DEz*Ex - DEx*Ez;
% C_Ez=DEx*Ey - DEy*Ex;
% 
% C_Hx=DHy*Hz - DHz*Hy;
% C_Hy=DHz*Hx - DHx*Hz;
% C_Hz=DHx*Hy - DHy*Hx;
% 
% Hx=-c0/mu_xx*C_Ex;
% Hy=-c0/mu_yy*C_Ey;
% Hz=-c0/mu_zz*C_Ez;
% 
% Ex=c0/eps_xx*C_Hx;
% Ey=c0/eps_yy*C_Hy;
% Ez=c0/eps_zz*C_Hz;

clear all
close all
Lz=1;
N=2^7
hz=Lz/N;
tau=1*hz;
tmax=10000;
c0=1;
epsvac=1;
letmebePeriodic=1;
%Hx, Ey

%ICS;
Hx=zeros([tmax N]);
Curl_Hx=zeros([tmax N]);
Ey=zeros([tmax N]);
Curl_Ey=zeros([tmax N]);

geom=@(x) (x>N/3).*(x<2*N/3);

mu=ones([1,N]);
eps=ones([1,N]);%(1-geom(1:N))+2*geom(1:N);

m_Hx=c0*tau./mu;
m_Ey=c0*tau./eps;

fig2=figure(2)
ax2=axes();
axis tight
set(ax2,'nextplot',...
'replacechildren');


aux_m_diag=diag(ones([1, N]),0)
aux_m_diag2=diag(ones([1, N]),1);
aux_m_diag2=aux_m_diag2(1:N,1:N);
aux_m_diag3=diag(ones([1, N]),-1);
aux_m_diag3=aux_m_diag3(1:N,1:N);

Dforw=aux_m_diag2-aux_m_diag
Dback=-aux_m_diag3+aux_m_diag

Dback(1,N)=-letmebePeriodic;
Dforw(N,1)=letmebePeriodic;

% weirdly enough, time should be in the first coordinate.
% otherwise everything is slower
for t=2:tmax
    t
    plot(Ey);
    
    shading interp
    %Hx(z, t+1/2tau)
    Curl_Ey=Dforw*( Ey(t-1,:)')/hz;
    Hx(t,:)=Hx(t-1,:)+(diag(m_Hx)*Curl_Ey)';
    Curl_Hx=Dback*( Hx(t  ,:)')/hz;
    Ey(t,:)=Ey(t-1,:)+(diag(m_Ey)*Curl_Hx)';
    
    %!!! it is VERY important to use Hx(t) from THIS time step
    %diverges otherwise
    
    Ey(t,floor(N/3))=Ey(t,floor(N/3))+gaussian_source(t,20,3);
    
figure(2);

plot(Ey(t,:));
hold on
plot(Hx(t,:));
hold off
axis([1 N -1.2 1.2])
drawnow;
end