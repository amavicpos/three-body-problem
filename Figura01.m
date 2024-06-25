%% Figure 01: Elliptic orbits, real Kepler-16 system
%IS units unless otherwise stated
tic
clear all ; clc ;

%% Kepler-16 data (wiki)
AU=1.496E11 ; MS=1.988E30 ; MJ=1.898E27 ; day=86400 ;  % units, AU: Astronomica Unit
m1=0.6897*MS ; m2=0.20255*MS ; d=0.22*AU ; T=41.079*day ; % orbit
r1=d*m2/(m1+m2) ; r2=d*m1/(m1+m2) ; % get stars positions
v1=2*pi*r1/T ; v2=2*pi*r2/T ; % get stars velocities
m3=0.333*MJ ; r3=0.7048*AU ; T3=228.776*day ; v3=2*pi*r3/T3 ; % planet 

%% DEFINE BODIES
m=[m1 m2 m3]; % masses
r0=[[-r1 0 0]' [+r2 0 0]' [r3 0 0]']; % initial r
v0=[[0 -v1 0]' [0 +v2 0]' [0 +v3 0]']; % initial v

%%
tend=T3*2 ; % final time
NS=5E5 ; % number of steps
[vx,vy,vz,x,y,z,K,U,Ktot,Utot,E,Ltot,Lxtot,Lytot,Lztot,t,L1,L2,L3]=VerletGravity(m,r0,v0,r3,tend,NS);

%% PLOTS
subplot(1,2,1) ; plot(x(1,:)/AU,y(1,:)/AU,x(2,:)/AU,y(2,:)/AU,x(3,:)/AU,y(3,:)/AU) ;
xlabel('Posición x (AU)'); ylabel('Posición y (AU)');
title('Trayectoria del sistema Kepler-16'); axis equal ; axis square ;
legend('Estrella Kepler-16A','Estrella Kepler-16B','Planeta Kepler-16b');
subplot(1,2,2) ; plot(t/day,Ktot,t/day,Utot,t/day,E); xlabel('Tiempo (días)');
ylabel('Energía (J)'); legend('Ktot','Utot','E'); title('Energías del sistema');
toc