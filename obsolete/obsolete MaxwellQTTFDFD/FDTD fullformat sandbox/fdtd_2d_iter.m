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
CEy=zeros([tmax N]);

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
    for ix=1:Nx   
        for iy=1:Ny
            CEx(ix,iy) = (Ez(t-1,ix,iy+1) - Ez(t-1,ix,iy))/hy;
            CEy(ix,iy) = -(Ez(t-1,ix+1,iy) - Ez(t-1,ix,iy))/hx; 
            Hx(t,ix,iy)=Hx(t-1,iz)+m_Hx(iz)*CEx(ix,iy);
            Hy(t,ix,iy)=Hy(t-1,iz)+m_Hy(iz)*CEy(ix,iy);
        end
    end

    
    for iy=1:Ny
        for ix=1:Nx
            CHz(ix,iy) = (Hy(t,ix,iy) - Hy(t,ix-1,iy))/hx-...
                (Hx(t,ix,iy) - Hx(t,ix,iy-1))/hy;
            Ez(t,ix,iy)=Ez(t-1,ix,iy)+m_Hy(iz)*CHz(ix,iy);
        end
    end
figure(2);
plot(Ey(t,:));
hold on
plot(Hx(t,:));
hold off
axis([1 N -1.2 1.2])
drawnow;
end