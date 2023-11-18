
pkg load control

g = (393/1025)*9.8;
m=360*(1025/170);

Cd=1.05;
rho = 0.015;
A=1.56*1.5*2+2.7*3;
beta = Cd*rho*A;


#FpI = RFp;

R= Euler2R([1*pi/180, -2*pi/180, 0]);

ang =20*(pi/180);

pkg load symbolic
syms T1 T2 T3 T4;

T11(T1) = [0; -sin(ang); -cos(ang)]*T1; #os rockets estam fixos, não temos esse estado relativo aos angulos dos rockets
T22(T2) = [0; sin(ang); -cos(ang)]*T2;
T33(T3) = [0; sin(ang); -cos(ang)]*T3;
T44(T4) = [0; -sin(ang); -cos(ang)]*T4;

FpI(T1,T2,T3,T4)= R*(T11(T1)+T22(T2)+T33(T3)+T44(T4));
Fg = m*g*[0; 0; 1];

Tum=20000;
Tdois(T3,T4) = solve((FpI(Tum,T2,T3,T4)-Fg)(3) ==[0;0;0](3));
Ttres(T4) = solve((FpI(Tum,Tdois(T3,T4),T3,T4)-Fg)(2) ==[0;0;0](2));
Tquatro = solve((FpI(Tum,Tdois(Ttres(T4),T4),Ttres(T4),T4)-Fg)(1) ==[0;0;0](1));
round(Tquatro)
round(Ttres(Tquatro))
round(Tdois(Ttres(Tquatro),Tquatro))
Tum


#FpI(2,2,3,10)



