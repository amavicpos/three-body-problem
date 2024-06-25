%% Figure 06 with CM2
%IS units unless otherwise stated
tic
clear all ; clc ; G=6.674E-11 ;

%% Kepler-16 data (wiki)
AU=1.496E11 ; MS=1.988E30 ; MJ=1.898E27 ; day=86400 ;  % units, AU: Astronomica Unit
m1=0.6897*MS ; m2=m1 ; T=2*41.079*day ; d=0.22*2*AU; % orbit % !!! Different T, d
r1=d*m2/(m1+m2) ; r2=d*m1/(m1+m2) ; % get stars positions
v1=2*pi*r1/T ; v2=2*pi*r2/T ; % get stars velocities
m3=0.63*MJ ; r3=0.22*2*AU ; T3=T*0.42 ; v3=2*pi*r3/T3 ; % planet % !!! Different m3, r3

%% DEFINE BODIES
m=[m1 m2 m3]; % masses
r0=[[0.25*r1  -0.01*r1 0]' [-0.25*r2 0*r2 0]' [-0.15*r3 0*r3 0]']; % initial r % !!! Different r0, v0
v0=[[-0.02*v1 1.37*v1 0]' [-0.1*v2 -2.12*v2 0]' [0.12*v3 0.75*v3 0]']; % initial v

%% FIND static and zero Center Of Mass 
N=length(m) ; % find number of bodies
rCM=sum(m.*r0,2)/sum(m) ;
vCM=sum(m.*v0,2)/sum(m) ;
r00=r0-rCM ; v00=v0-vCM ; % CORRECT POSITIONS AND VELOCITIES
plot(r00(1,:),r00(2,:),'o') ; 
rmax=r3 ; xlim([-rmax,rmax]) ; ylim([-rmax,rmax]) ; axis square ; grid on ;

%% SET INTEGRATION PARAMETERS
tend=T3*2 ; % final time % !!! Different tend
NS=5E5 ; % number of steps
tfinal=round(tend); dt=tfinal/NS ; t=1:dt:tfinal ;

%% INTEGRATE
Gmm=-G*m.*m' ; Gmm2=-G*m(1:2:3).*m(1:2:3)' ; %Gmm2 for stars' Center of Mass
r000=[norm(r00(:,1)-r00(:,2)) norm(r00(:,1)-r00(:,3)) norm(r00(:,2)-r00(:,3))]';
r(:,1)=r000; %relative initial positions
K00=0.5*m.*[norm(v00(:,1)) norm(v00(:,2)) norm(v00(:,3))].^2 ;
U00(1)=Gmm(1,2)/r000(1)+Gmm(1,3)/r000(2);
U00(2)=Gmm(2,1)/r000(1)+Gmm(2,3)/r000(3);
U00(3)=Gmm(3,1)/r000(2)+Gmm(3,2)/r000(3);
[x,y,z,vx,vy,vz]=deal(zeros(N,NS)) ;
[x2,y2,z2,vx2,vy2,vz2]=deal(zeros(2,NS)) ;
[K,U,L1,L2,L3]=deal(zeros(3,NS)) ; [Ktot,Utot,E,Ltot]=deal(zeros(1,NS)) ;
x(:,1)=r00(1,:)' ; y(:,1)=r00(2,:)' ; z(:,1)=r00(3,:)' ; %separate coordinates
vx(:,1)=v00(1,:)' ; vy(:,1)=v00(2,:)' ; vz(:,1)=v00(3,:)' ;
x2(:,1)=r00(1,2:3)' ; y2(:,1)=r00(2,2:3)' ; z2(:,1)=r00(3,2:3)' ;
vx2(:,1)=v00(1,2:3)' ; vy2(:,1)=v00(2,2:3)' ; vz2(:,1)=v00(3,2:3)' ;
K(:,1)=K00' ; U(:,1)=U00' ;
L1(:,1)=cross([x(1,1) y(1,1) z(1,1)],m(1).*[vx(1,1) vy(1,1) vz(1,1)])' ;
L2(:,1)=cross([x(2,1) y(2,1) z(2,1)],m(2).*[vx(2,1) vy(2,1) vz(2,1)])';
L3(:,1)=cross([x(3,1) y(3,1) z(3,1)],m(3).*[vx(3,1) vy(3,1) vz(3,1)])';
Ltot(1)=norm(L1(:,1)+L2(:,1)+L3(:,1));

for i=1:NS-1
   ddx=x(:,i)-x(:,i)' ; ddy=y(:,i)-y(:,i)' ; ddz=z(:,i)-z(:,i)' ;
   dd=sqrt(ddx.^2+ddy.^2+ddz.^2) ;
   ffx=Gmm.*ddx./dd.^3 ; ffy=Gmm.*ddy./dd.^3 ; ffz=Gmm.*ddz./dd.^3 ;
   ffx(1:N+1:end)=0 ; ffy(1:N+1:end)=0 ; ffz(1:N+1:end)=0 ; %ignore force over oneself
   fx=sum(ffx,2) ; fy=sum(ffy,2) ; fz=sum(ffz,2) ; %sum rows, force over each body
   
   vx(:,i+1)=vx(:,i)+(fx./(m'))*dt ;
   vy(:,i+1)=vy(:,i)+(fy./(m'))*dt ;
   vz(:,i+1)=vz(:,i)+(fz./(m'))*dt ;
   x(:,i+1)=x(:,i)+vx(:,i+1)*dt ;
   y(:,i+1)=y(:,i)+vy(:,i+1)*dt ;
   z(:,i+1)=z(:,i)+vz(:,i+1)*dt ;
   v1=norm([vx(1,i+1) vy(1,i+1) vz(1,i+1)]);
   v2=norm([vx(2,i+1) vy(2,i+1) vz(2,i+1)]);
   v3=norm([vx(3,i+1) vy(3,i+1) vz(3,i+1)]);
   vx2(:,i+1)=vx2(:,i)+(fx(2:3)./(m(2:3)'))*dt ; %Two stars' vectors
   vy2(:,i+1)=vy2(:,i)+(fy(2:3)./(m(2:3)'))*dt ;
   vz2(:,i+1)=vz2(:,i)+(fz(2:3)./(m(2:3)'))*dt ;
   x2(:,i+1)=x2(:,i)+vx2(:,i+1)*dt ;
   y2(:,i+1)=y2(:,i)+vy2(:,i+1)*dt ;
   z2(:,i+1)=z2(:,i)+vz2(:,i+1)*dt ;
   rCM2(:,i+1)=sum(m(2:3).*[x2(:,i+1)'; y2(:,i+1)'; z2(:,i+1)'],2)/sum(m(2:3)) ;
   vCM2(:,i+1)=sum(m(2:3).*[vx2(:,i+1)'; vy2(:,i+1)'; vz2(:,i+1)'],2)/sum(m(2:3)) ;

   K(:,i+1)=(0.5*m.*[v1 v2 v3].^2)' ;
   r(1,i+1)=norm([x(1,i+1)-x(2,i+1) y(1,i+1)-y(2,i+1) z(1,i+1)-z(2,i+1)]);
   r(2,i+1)=norm([x(1,i+1)-x(3,i+1) y(1,i+1)-y(3,i+1) z(1,i+1)-z(3,i+1)]);
   r(3,i+1)=norm([x(2,i+1)-x(3,i+1) y(2,i+1)-y(3,i+1) z(2,i+1)-z(3,i+1)]);
   U(1,i+1)=Gmm(1,2)/r(1,i+1)+Gmm(1,3)/r(2,i+1);
   U(2,i+1)=Gmm(2,1)/r(1,i+1)+Gmm(2,3)/r(3,i+1);
   U(3,i+1)=Gmm(3,1)/r(2,i+1)+Gmm(3,2)/r(3,i+1);
   L1(:,i+1)=cross([x(1,i+1) y(1,i+1) z(1,i+1)],m(1).*[vx(1,i+1) vy(1,i+1) vz(1,i+1)])';
   L2(:,i+1)=cross([x(2,i+1) y(2,i+1) z(2,i+1)],m(2).*[vx(2,i+1) vy(2,i+1) vz(2,i+1)])';
   L3(:,i+1)=cross([x(3,i+1) y(3,i+1) z(3,i+1)],m(3).*[vx(3,i+1) vy(3,i+1) vz(3,i+1)])';
   Ltot(i+1)=norm(L1(:,i+1)+L2(:,i+1)+L3(:,i+1));
   
   if rem(i,5000)==0 %Animated plot
       subplot(1,2,1) ;
       plot(x(:,i)/AU,y(:,i)/AU,'o',0,0,'xm',rCM2(1,i+1)/AU,rCM2(2,i+1)/AU,'xg') ;
       title('Evolución del sistema'); legend('Estrellas', 'Centro de masas total', 'CM estrellas entrelazadas');
       rmax=1.5*r3/AU ; xlim([-rmax,rmax]) ; ylim([-rmax,rmax]) ;
       xlabel('Posición x (AU)') ; ylabel('Posición y (AU)') ; axis square ; grid on ;
       drawnow ;
   end    
end

Lxtot=L1(1,:)+L2(1,:)+L3(1,:);
Lytot=L1(2,:)+L2(2,:)+L3(2,:);
Lztot=L1(3,:)+L2(3,:)+L3(3,:);

%% PLOTS
Ktot=sum(K,1) ; Utot=0.5*sum(U,1); E=Ktot+Utot ;
subplot(1,3,1) ; plot(x(1,:)/AU,y(1,:)/AU,x(2,:)/AU,y(2,:)/AU,x(3,:)/AU,y(3,:)/AU) ;
hold on; plot(0,0,'xm',rCM2(1,end)/AU,rCM2(2,end)/AU,'xg'); hold off;
xlabel('Posición x (AU)'); ylabel('Posición y (AU)');
title('Trayectoria del Sistema Kepler-16 modificado'); axis equal ; axis square ;
legend('Estrella Kepler-16A','Estrella Kepler-16A','Estrella Kepler-16A','Centro de masas total','CM estrellas entrelazadas');
subplot(1,3,2) ; plot(t/day,Ktot,t/day,Utot,t/day,E); xlabel('Tiempo (días)');
ylabel('Energía (J)'); legend('Ktot','Utot','E'); title('Energías del Sistema');
subplot(1,3,3); plot(t/day,Ltot,t/day,Lxtot,t/day,Lytot,t/day,Lztot);
title('Momento angular total'); xlabel('Tiempo (días)');
ylabel('Momento angular ((kg·m^2)/s)'); legend('Ltot (módulo)','Lxtot','Lytot','Lztot');
toc