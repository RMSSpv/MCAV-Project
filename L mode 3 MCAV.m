clear all;


%initPlots;

g = (393/1025)*9.8;
mL=360*(1025/170); %360kg é a masssa do Lander de outra missão e 170kg do Rover, decidi manter a proporção
mR=1025;
m=mR+mL;


Cd=1.05;
rho = 0.015;
A=diag([13.08 , 7.62 , 12.84]);
beta = Cd*rho*A;


CML=[1.35; 3; -2.74];
CMR=[1.35; 3; -1.1];
CM = (mL*CML+mR*CMR)*(1/m);
CMLP1 = CML +[0.39; 1.5; 0];
CMLP2 = CML +[-0.39; -1.5; 0];
CMLP3 = CML +[0.39; -1.5; 0];
CMLP4 = CML +[-0.39; 1.5; 0];
CMRP1 = CMR +[0.675; 0.75; 0];
CMRP2 = CMR +[-0.675; -0.75; 0];
CMRP3 = CMR +[0.675; -0.75; 0];
CMRP4 = CMR +[-0.675; 0.75; 0];
JL1 = ((mL/4)*((CMLP1-CM).^2));
JL2 = ((mL/4)*((CMLP2-CM).^2));
JL3 = ((mL/4)*((CMLP3-CM).^2));
JL4 = ((mL/4)*((CMLP4-CM).^2));
JR1 = ((mR/4)*((CMRP1-CM).^2));
JR2 = ((mR/4)*((CMRP2-CM).^2));
JR3 = ((mR/4)*((CMRP3-CM).^2));
JR4 = ((mR/4)*((CMRP4-CM).^2));
Jp=(JL1+JR1+JL2+JR2+JL3+JR3+JL4+JR4);
Jx = Jp(1) ;
Jy=Jp(2);
Jz=Jp(3);
J = diag([ Jx , Jy , Jz ]);


Dt = 0.01;
t = 0:Dt:50;
N = length( t );

T1 = 12777.57858; %T de eq e inicial
T2 = 12777.57858;
T3 = 12777.57858;
T4 = 12777.57858;
ang = 20*pi/180 ;
T1u = (t>=0)*T1;
T2u = (t>=0)*T2;
T3u = (t>=0)*T3;
T4u = (t>=0)*T4;
u=[T1u;T2u;T3u;T4u];

r1=[(2.7+1.56)/2;6;-2.2]-CM;
r2=[(2.7-1.56)/2;0;-2.2]-CM;
r3=[(2.7+1.56)/2;0;-2.2]-CM;
r4=[(2.7-1.56)/2;6;-2.2]-CM;

r1x = r1(1);
r1y = r1(2);
r1z = r1(3);
r2x = r2(1);
r2y = r2(2);
r2z = r2(3);
r3x = r3(1);
r3y = r3(2);
r3z = r3(3);
r4x = r4(1);
r4y = r4(2);
r4z = r4(3);


vx = 0 ;
vy = 0 ;
vz = 0 ;
fi = 0 ;
theta = 0 ;
psi = 0 ;
px = 0 ;
py = 0 ;
pz = - 21.3 ;
wx =0 ;
wy =0 ;
wz = 0 ;
x0 = [px ; py ; pz ; fi; theta ; psi ; vx ; vy ; vz ; wx ; wy ; wz ];

A1 = zeros (3,3) ; %f1/p

A2 = [ vy*(cos (fi)*sin(theta)*cos(psi)+sin(fi)*sin(psi))+vz*(-sin(fi)*sin(theta)*cos(psi)+cos(fi)*sin(psi)) , vx*(-sin(theta)*cos(psi))+vy*(sin(fi)*cos(theta)*cos(psi))+vz*(cos(fi)*cos(theta)*cos(psi)) , vx*(-cos(theta)*sin(psi))+vy*(-sin(fi)*sin(theta)*sin(psi)-cos(fi)*cos(psi))+vz*(-cos(fi)*sin(theta)*sin(psi)+sin(fi)*cos(psi)) ;
       vy*(cos(fi)*sin(theta)*sin(psi)-sin(fi)-cos(psi))+vz*(-sin(fi)*sin(theta)*sin(psi)-cos(fi)*cos(psi)) , vx*(-sin(theta)*sin(psi))+vy*(sin(fi)*cos(theta)*sin(psi))+vz*(cos(fi)*cos(theta)*sin(psi)) , vx*(cos(theta)*cos(psi))+vy*(sin(fi)*sin(theta)*cos(psi)-cos(fi)*sin(psi))+vz*(cos(fi)*sin(theta)*cos(psi)+sin(fi)*sin(psi)) ;
       cos(theta)*(vy*cos(fi)-vz*sin(fi)) , -vx*cos(theta)+vy*(sin(fi)*sin(theta))+vz*(cos(fi)*sin(theta)) , 0 ] ; %f1/lambda

A3 = [ cos(theta)*cos(psi) , sin(fi)*sin(theta)*cos(psi)-cos(fi)*sin(psi) , cos(fi)*sin(theta)*cos(psi)+sin(fi)*sin(psi) ;
       cos(theta)*sin(psi) , sin(fi)*sin(theta)*sin(psi)+cos(fi)*cos(psi) , cos(fi)*sin(theta)*sin(psi)-sin(fi)*cos(psi) ;
       -sin(theta) , sin(fi)*cos(theta) , cos(fi)*cos(theta) ] ; %f1/v

A4 = zeros (3,3) ; %f1/w

A5 =  zeros (3,3) ; %f2/p

A6 = [wy*tan(theta)*cos(fi) - wz*tan(theta)*sin(fi) , wy*sin(fi)* (sec (theta))^2 + wz*cos(fi)*(sec (theta))^2 , 0 ;
      -wy*(sin(fi)-wz*cos(theta)) , 0 , 0 ;
      1/cos(theta)*(wy*cos(fi)-wz*sin(fi)) , (sec(theta)*tan(theta))*(wy*sin(fi)+wz*cos(theta)) , 0 ] ; %f2/lambda

A7 = zeros (3,3) ; %f2/v

A8 = [ 1 , sin(fi)*tan(theta) , cos(fi)*tan(theta) ;
       0 , cos(fi) , -sin(fi) ;
       0 , sin(fi)/cos(theta) , cos(fi)/cos(theta) ] ; % f2/w

A9 = zeros (3,3) ; %f3/p

A10 = [0 , -cos(theta)*g, 0;
       cos(fi)*cos(theta)*g, -sin(theta)*sin(fi)*g, 0;
       -sin(fi)*cos(theta)*g, -sin(theta)*cos(fi)*g, 0]; %f3/lambda

A11 = [2*beta(1,1)*vx*(1/m), wz, -wy;
       -wz, 2*beta(2,2)*vy*(1/m), wx;
       wy, -wx, 2*beta(3,3)*vz*(1/m)]; %f3/v

A12 = [0, -vz, vy;
       vz, 0, -vx;
       -vy, vx, 0]; %f3/w


A13 = zeros(3,3);

A14 = zeros(3,3);

A15 = zeros(3,3);

A16 = [0, wz*(Jy-Jz)*(1/Jx), wy*(Jy-Jz)*(1/Jx);
wz*(Jz-Jx)*(1/Jy), 0, wx*(Jz-Jx)*(1/Jy);
wy*(Jx-Jy)*(1/Jz), wx*(Jx-Jy)*(1/Jz), 0] ;

A = [A1 , A2 , A3 , A4 ;
     A5 , A6 , A7 , A8 ;
     A9 , A10 , A11 , A12 ;
     A13 , A14 , A15 , A16 ] 
 
 
B1 = zeros (3,1) ;

B2 = zeros (3,1) ;

B3 = zeros (3,1) ;

B4 = zeros (3,1) ;

B5 =  zeros (3,1) ;

B6 =  zeros (3,1) ;

B7 =  zeros (3,1) ;

B8 =  zeros (3,1) ;

B9 = 1/m*[0 ; -T1*sin(ang) ; -T1*cos(ang)] ;

B10 = 1/m*[0 ; T2*sin(ang) ; -T2*cos(ang)] ;

B11 = 1/m*[0 ; T3*sin(ang) ; -T3*cos(ang)] ;

B12 =  1/m*[0 ; -T4*sin(ang) ; -T4*cos(ang)] ;

B13 = [-T1/Jx*(r1z*sin(ang)*(-1)+r1y*cos(ang)) ; T1/Jy * (r1x*cos(ang)) ; T2/Jz * (r2x*sin(ang)*(-1)) ] ;

B14 = [-T2/Jx*(r2z*sin(ang)+r2y*cos(ang)) ; T2/Jy * (r2x*cos(ang)) ; T2/Jz * (r2x*sin(ang)) ] ;

B15 = [-T3/Jx*(r3z*sin(ang)*+r3y*cos(ang)) ; T3/Jy * (r1x*cos(ang)) ; T3/Jz * (r3x*sin(ang)) ] ;

B16 =[-T4/Jx*(r4z*sin(ang)*(-1)+r4y*cos(ang)) ; T4/Jy * (r4x*cos(ang)) ; T4/Jz * (r4x*sin(ang)*(-1)) ] ;

B = [B1 , B2 , B3 , B4 ;
     B5 , B6 , B7 , B8 ;
     B9 , B10 , B11 , B12 ;
     B13 , B14 , B15 , B16 ] 

C = eye(12);

D = zeros(12,4);

sys = ss(A,B,C,D);


y_L = lsim(sys,u,t,x0)';

[Vj,Jor] = jordan(A)
[V,DL,W] = eig(A);
mode_obs = C*V
mode_ctrl = W'*B


