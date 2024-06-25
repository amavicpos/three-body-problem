%% Extra figure Verlet: Three bodies with the same mass Verlet
%IS units unless otherwise stated
tic
clear all ; clc ; G=6.674E-11 ;

%% Kepler-16 data (wiki)
AU=1.496E11 ; MS=1.988E30 ; MJ=1.898E27 ; day=86400 ;  % units, AU: Astronomica Unit
m1=0.6897*MS ; m2=m1 ; d=0.22*2*AU ; m3=m1 ; T=41.079*day ; % orbit % !!! Different m2, d, m3
r1=d*m2/(m1+m2) ; r2=d*m1/(m1+m2) ; % get stars positions
v1=2*pi*r1/T ; v2=2*pi*r2/T ; % get stars velocities
r3=0.7048*AU ; T3=228.776*day ; v3=2*pi*r3/T3*0.9 ; % planet % !!! Different v3

%% DEFINE BODIES
m=[m1 m2 m3]; % masses
r0=[[0.2693*r3 -1.0020*r3 0]' [-0.2328*r3 -0.5978*r3 0]' [-0.0347*r3 1.1856*r3 0]']; % initial r % !!! Different r0, v0
v0=[[0.2059*v3 -0.9396*v3 0]' [-0.4553*v3 1.0471*v3 0]' [0.2495*v3 -0.1076*v3 0]']; % initial v

%%
tend=T3*130 ; % final time % !!! Different tend
NS=5E5 ; % number of steps
[vx,vy,vz,x,y,z,K,U,Ktot,Utot,E,Ltot,Lxtot,Lytot,Lztot,t]=VerletGravity(m,r0,v0,r3,tend,NS);

%% PLOTS
Ktot=sum(K,1) ; Utot=0.5*sum(U,1); E=Ktot+Utot ;
subplot(1,3,1) ; plot(x(1,:)/AU,y(1,:)/AU,x(2,:)/AU,y(2,:)/AU,x(3,:)/AU,y(3,:)/AU) ;
xlabel('Posición x (AU)'); ylabel('Posición y (AU)');
title('Trayectoria del Sistema Kepler-16 modificado'); axis equal ; axis square ;
legend('Estrella Kepler-16A','Estrella Kepler-16A','Estrella Kepler-16A');
subplot(1,3,2) ; plot(t,Ktot,t,Utot,t,E); xlabel('Tiempo (s)');
ylabel('Energía (J)'); legend('Ktot','Utot','E'); title('Energías del Sistema');
subplot(1,3,3); plot(t,Ltot,t,Lxtot,t,Lytot,t,Lztot);
title('Momento angular total, módulo'); xlabel('Tiempo (s)');
ylabel('Momento angular ((kg·m^2)/s)'); legend('Ltot (módulo)','Lxtot','Lytot','Lztot');
toc