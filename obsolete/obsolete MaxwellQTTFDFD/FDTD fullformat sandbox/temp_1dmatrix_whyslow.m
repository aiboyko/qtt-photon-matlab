clear all
close all
Lz=1;
Nz=2^7
hz=Lz/Nz;
tau=1*hz;
tmax=1000;
c0=1;
epsvac=1;
letmebePeriodic=1;
%Hx, Ey

%ICS;
Hx=zeros([ Nz tmax]);
C_Hx=zeros([ Nz tmax]);
Ey=zeros([ Nz tmax]);
C_Ey=zeros([ Nz tmax]);

geom=@(x) (x>Nz/3).*(x<2*Nz/3);

mu=ones([1,Nz]);
eps=(1-geom(1:Nz))+3*geom(1:Nz);

m_Hx=c0*tau./mu;
m_Ey=c0*tau./eps;

fig2=figure(2)
ax2=axes();
axis tight
set(ax2,'nextplot',...
'replacechildren');


aux_m_diag=diag(ones([1, Nz]),0)
aux_m_diag2=diag(ones([1, Nz]),1);
aux_m_diag2=aux_m_diag2(1:Nz,1:Nz);
aux_m_diag3=diag(ones([1, Nz]),-1);
aux_m_diag3=aux_m_diag3(1:Nz,1:Nz);

Dforw=aux_m_diag2-aux_m_diag
Dback=-aux_m_diag3+aux_m_diag

Dback(1,Nz)=-letmebePeriodic;
Dforw(Nz,1)=letmebePeriodic;


for t=2:tmax
    t
    plot(Ey);
    
    
    %Hx(z, t+1/2tau)
    
    C_Ex=Dforw*(Ey(:,t-1))/hz;
    Hx(:,t)=Hx(:,t-1)+diag(m_Hx)*C_Ex;
    
    C_Hx=Dback*(Hx(:,t))/hz;
    Ey(:,t)=Ey(:,t-1)+diag(m_Ey)*C_Hx;
    %!!! it is VERY important to use Hx(t) from THIS time step
    %diverges otherwise
    
    Ey(floor(Nz/3),t)=Ey(floor(Nz/3),t)+gaussian_source(t,20,3);
    
figure(2);

plot(Ey(:,t));
hold on
plot(Hx(:,t));
hold off
axis([1 Nz -1.2 1.2])
drawnow;
end