%% Figure 07: Noncoplanar orbits
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
r0=[[-1.3865e+09 -2.8330e+09 -4.8565e+09]' [4.8948e+09 9.5528e+09 1.6524e+10]' [-1.1054e+11 5.9780e+10 8.1163e+09]']; % initial r
v0=[[2.3869e+03 -1.1027e+04 8.1174e+03]' [-8.1045e+03 3.7587e+04 -2.7640e+04]' [-1.4722e+04 -2.5759e+04 -219.5053]']; % initial v

%%
tend=T3*1.3 ; % final time
NS=5E5 ; % number of steps
[vx,vy,vz,x,y,z,K,U,Ktot,Utot,E,Ltot,Lxtot,Lytot,Lztot,t]=VerletGravity(m,r0,v0,r3,tend,NS);

%% PLOTS
subplot(2,1,1) ;
plot3(x(1,:)/AU,y(1,:)/AU,z(1,:)/AU,x(2,:)/AU,y(2,:)/AU,z(2,:)/AU,x(3,:)/AU,y(3,:)/AU,z(3,:)/AU) ;
xlabel('Posición x (AU)'); ylabel('Posición y (AU)'); zlabel('Posición z (AU)');
title('Trayectoria del sistema Kepler-16 modificado');
legend('Estrella Kepler-16A','Estrella Kepler-16B','Planeta Kepler-16b');
subplot(2,1,2) ; plot(t/day,Ktot,t/day,Utot,t/day,E); xlabel('Tiempo (días)');
ylabel('Energía (J)'); legend('Ktot','Utot','E'); title('Energías del sistema');
toc