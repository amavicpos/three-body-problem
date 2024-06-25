function [vx,vy,vz,x,y,z,K,U,Ktot,Utot,E,Ltot,Lxtot,Lytot,Lztot,t,L1,L2,L3]=EulerGravity(m,r0,v0,r3,tend,NS)
G=6.674E-11 ; AU=1.496E11 ;
%% FIND static and zero Center Of Mass 
N=length(m) ; % find number of bodies
rCM=sum(m.*r0,2)/sum(m) ;
vCM=sum(m.*v0,2)/sum(m) ;
r00=r0-rCM ; v00=v0-vCM ; % CORRECT POSITIONS AND VELOCITIES

%% SET INTEGRATION PARAMETERS
tfinal=round(tend); dt=tfinal/NS ; t=1:dt:tfinal ;

%% INTEGRATE
Gmm=-G*m.*m' ;
r000=[norm(r00(:,1)-r00(:,2)) norm(r00(:,1)-r00(:,3)) norm(r00(:,2)-r00(:,3))]';
r(:,1)=r000; %relative initial positions
K00=0.5*m.*[norm(v00(:,1)) norm(v00(:,2)) norm(v00(:,3))].^2 ;
U00(1)=Gmm(1,2)/r000(1)+Gmm(1,3)/r000(2);
U00(2)=Gmm(2,1)/r000(1)+Gmm(2,3)/r000(3);
U00(3)=Gmm(3,1)/r000(2)+Gmm(3,2)/r000(3);
[x,y,z,vx,vy,vz]=deal(zeros(N,NS)) ;
[K,U,L1,L2,L3]=deal(zeros(3,NS)) ; [Ktot,Utot,E,Ltot]=deal(zeros(1,NS)) ;
x(:,1)=r00(1,:)' ; y(:,1)=r00(2,:)' ; z(:,1)=r00(3,:)' ; %separate coordinates
vx(:,1)=v00(1,:)' ; vy(:,1)=v00(2,:)' ; vz(:,1)=v00(3,:)' ;
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
       plot(x(:,i)/AU,y(:,i)/AU,'o') ; title('Evolución del sistema');
       rmax=1.5*r3/AU ; xlim([-rmax,rmax]) ; ylim([-rmax,rmax]) ;
       xlabel('Posición x (AU)') ; ylabel('Posición y (AU)') ; axis square ; grid on ;
       drawnow ;
   end    
end

Lxtot=L1(1,:)+L2(1,:)+L3(1,:);
Lytot=L1(2,:)+L2(2,:)+L3(2,:);
Lztot=L1(3,:)+L2(3,:)+L3(3,:);
Ktot=sum(K,1) ; Utot=0.5*sum(U,1); E=Ktot+Utot ;
end