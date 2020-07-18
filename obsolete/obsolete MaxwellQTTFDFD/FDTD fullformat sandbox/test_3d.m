
c0=1;


C_Ex=DEy*Ez - DEz*Ey;
C_Ey=DEz*Ex - DEx*Ez;
C_Ez=DEx*Ey - DEy*Ex;

C_Hx=DHy*Hz - DHz*Hy;
C_Hy=DHz*Hx - DHx*Hz;
C_Hz=DHx*Hy - DHy*Hx;

Hx=-c0/mu_xx*C_Ex;
Hy=-c0/mu_yy*C_Ey;
Hz=-c0/mu_zz*C_Ez;

Ex=c0/eps_xx*C_Hx;
Ey=c0/eps_yy*C_Hy;
Ez=c0/eps_zz*C_Hz;