clear all
close all
Lz=1;
N=2^9
hz=Lz/N;
tau=1*hz;
tmax=10000;
c0=1;
epsvac=1;

ksrc=floor(N/5);

letmebePeriodic=1;

%Hx, Ey

%ICS;
Hx=zeros([tmax N]);
Curl_Hx=zeros([tmax N]);
Ey=zeros([tmax N]);
Curl_Ey=zeros([tmax N]);

geom=@(x) (x>N/3).*(x<2*N/3);

mu=ones([1,N]);
eps=ones([1,N]);%(1-geom(1:N))+3*geom(1:N);

m_Hx=c0*tau./mu;
m_Ey=c0*tau./eps;

fig2=figure(2)
ax2=axes();
axis tight
set(ax2,'nextplot',...
'replacechildren');

for t=2:tmax
    t   
    for iz=1:N
        if iz==N
            Curl_Ey(N) = (letmebePeriodic*Ey(t-1,1) - Ey(t-1,N))/hz;
        else
            Curl_Ey(iz) = (Ey(t-1,iz+1) - Ey(t-1,iz))/hz;
        end
        Hx(t,iz)=Hx(t-1,iz)+m_Hx(iz)*Curl_Ey(iz);
    end
    for iz=1:N
        if iz==1
            Curl_Hx(1) = (Hx(t,1) - letmebePeriodic*Hx(t,N))/hz;
        else
            Curl_Hx(iz) = (Hx(t,iz) - Hx(t,iz-1))/hz;
        end
        Ey(t,iz)=Ey(t-1,iz)+m_Ey(iz)*Curl_Hx(iz)+double(iz==ksrc)*gaussian_source(t,20,3);
    end

figure(2);
plot(Ey(t,:));
hold on
plot(Hx(t,:));
hold off
axis([1 N -1.2 1.2])
drawnow;
end