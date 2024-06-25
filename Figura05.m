%% Figure 05: Internal planetary orbit 01
%IS units unless otherwise stated
tic
clear all ; clc ;

%% Kepler-16 data (wiki)
AU=1.496E11 ; MS=1.988E30 ; MJ=1.898E27 ; day=86400 ;  % units, AU: Astronomica Unit
m1=0.6674*MS*1.3 ; m2=m1 ; d=0.22*2*AU ; T=2*41.079*day ; % orbit % !!! Different m2, d, T
r1=d*m2/(m1+m2) ; r2=d*m1/(m1+m2) ; % get stars positions
v1=2*pi*r1/T ; v2=2*pi*r2/T ; % get stars velocities
m3=0.6674*MJ*0.3 ; r3=0.64*0.22*2*AU ; T3=T*0.7 ; v3=2*pi*r3/T3 ; % planet % !!! Different m3, r3, T3

%% DEFINE BODIES
m=[m1 m2 m3]; % masses
r0=[[-3.2468E9  1.0689E10 0]' [3.2468E9 -1.0688E10 0]' [-5.9836E5 -2.4183E9 0]']; % initial r %Different r0, v0
v0=[[-4.8585E4 -1.4702E4 0]' [4.8549E4 1.4692E4 0]' [1.6529E5 4.4003E4 0]']; % initial v

%%
tend=T3*0.835 ; % final time
NS=5E5 ; % number of steps
[vx,vy,vz,x,y,z,K,U,Ktot,Utot,E,Ltot,Lxtot,Lytot,Lztot,t]=VerletGravity(m,r0,v0,r3,tend,NS);

%% PLOTS
subplot(1,2,1) ; plot(x(1,:)/AU,y(1,:)/AU,x(2,:)/AU,y(2,:)/AU,x(3,:)/AU,y(3,:)/AU) ;
xlabel('Posición x (AU)'); ylabel('Posición y (AU)');
title('Trayectoria del sistema Kepler-16 modificado'); axis equal ; axis square ;
legend('Estrella Kepler-16A','Estrella Kepler-16A','Planeta');
subplot(1,2,2) ; plot(t/day,Ktot,t/day,abs(Utot)); xlabel('Tiempo (días)');
ylabel('Energía (J)'); legend('Ktot','Utot (valor absoluto)'); title('Energías del sistema');
toc