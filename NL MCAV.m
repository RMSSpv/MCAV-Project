#f = @(x) x.^2 + 3*x - 1 + 5*x.*sin(x);
% these next lines take the Anonymous function into a symbolic formula
%pkg load symbolic
%syms x;
%ff = f(x);
% now calculate the derivative of the function
%ffd = diff(ff, x)


pkg load control
#graphics_toolkit('gnuplot')

%NOTA: o Eixo do Z do Inercial considero o virado para baixo

Dt = 0.01;

%Constantes que assumi
t = 0:Dt:50;
#g = 9.8;
g = (393/1025)*9.8

N = length( t );
#T=10; %estou a considerar um caso que o Thrust em todos os rotores é constante, e igual em todos os rotores(sem momento e conseuqnte variação de w ou x4)
%em uma nova versão quando tiver tempo meto a força de propulsão fp  e o momento np a influenciar o drone


x = zeros( 12 ,N+1);
u = zeros( 6 ,N);
y = zeros( 12 ,N);
fxu = zeros( 12 ,N);


mL=360*(1025/170); #360kg é a masssa do Lander de outra missão e 170kg do Rover, decidi manter a proporção
mR=1025;


m=mR+mL;
CML=[1.35; 3; -2.74];
CMR=[1.35; 3; -1.1];
CM = (mL*CML+mR*CMR)*(1/m); #Vamos considerar o centro de massa da estrutura toda o centro do Body

#posiçoes do Rockets relativamente ao centro de massa
r1=[(2.7+1.56)/2;6;-2.2]-CM;
r2=[(2.7-1.56)/2;0;-2.2]-CM;
r3=[(2.7+1.56)/2;0;-2.2]-CM;
r4=[(2.7-1.56)/2;6;-2.2]-CM;

#Para calcular o momento de Inercia vou dividir o lander em 4 paralelipedos iguais de massa igual a mL/4 e achar os respetivos centros de massa CMP em relação a CM, faço o mesmo para o Rover

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

Jx=(JL1+JR1+JL2+JR2+JL3+JR3+JL4+JR4)(1);
Jy=(JL1+JR1+JL2+JR2+JL3+JR3+JL4+JR4)(2);
Jz=(JL1+JR1+JL2+JR2+JL3+JR3+JL4+JR4)(3);

J = diag([Jx, Jy, Jz]);

#input u = [Fp + Fg; np]
ang = 20*(pi/180);

mod_T1 = ((1:N)>=0)*4800; #Rotor1 constante e igual a 5N para cada instante
mod_T2 = ((1:N)>=0)*4800; #Rotor2 constante e igual a 5N para cada instante
mod_T3 = ((1:N)>=0)*4800;
mod_T4 =((1:N)>=0)*4800;
u=[mod_T1; mod_T2; mod_T3; mod_T4]; #não meto os proprios T porque caso quisesse mudar o angulo dos rockets assim ja estaria bem e so seria adicionar os angulos aos inputs
T1 = [0; -sin(ang); -cos(ang)]*u(1,:); #os rockets estam fixos, não temos esse estado relativo aos angulos dos rockets
T2 = [0; sin(ang); -cos(ang)]*u(2,:);
T3 = [0; sin(ang); -cos(ang)]*u(3,:);
T4 = [0; -sin(ang); -cos(ang)]*u(4,:);

Fp = T1+T2+T3+T4; #Como é invariante no tempo pode estar fora do ciclo for se não tinha de estar la dentrro pois seria influenciado pelos outros estados, o mesmo serve para o np

ni1 = ni2 = ni3 = ni4 = [0; 0; 0];
Tor1 = skew(r1)*T1;
Tor2 = skew(r2)*T2;
Tor3 = skew(r3)*T3;
Tor4 = skew(r4)*T4;
Tor = Tor1 + Tor2 + Tor3+ Tor4;
ni = ni1 + ni2 + ni3+ ni4;
np = ni + Tor;


#u=[Fp; np]; #Devia ser Fp + Fg e acredito eu para considerar Fp + Fg, ou seja em vez de defenir  os valores dos rotores assumir que Fp + Fg vap tomar alguns valores



#Vamos fazer um modelo para quando este so desloca numa direção logo A = cte e Cd também

#Considerando o movimento descendente na direção, Z v = v_z
Cd=1.05;
rho = 0.015;
A=1.56*1.5*2+2.7*3;




%Os valores iniciais fui eu que assumi
x(1:3,1)=[0; 0; -2100];
x(4:6,1)=[0; 0; 0];
x(7:9,1)=[0; 0; 89];
x(10:12,1)=[0; 0; 0];

#R = Euler2R(x(4:6,1));
  #Forças
#Fg = m*g*R'*[0; 0; 1];
#aux = u(1:3,1)+Fg;

#mod_v=((x(7,1)^2)+(x(8,1)^2)+(x(9,1)^2))^(1/2);
#Fa = rho*(mod_v^2)*Cd*A*(1/2)*(-x(7:9,1)/mod_v);
#rrrr=(1/m)*Fa + (1/m)*aux;
for k = 1 :N

  R = Euler2R(x(4:6,k));
  #Forças
  Fg = m*g*R'*[0; 0; 1];
  #aux = u(1:3,k)+Fg; #??????
  #u(1:3,k) = aux;
  mod_v=((x(7,k)^2)+(x(8,k)^2)+(x(9,k)^2))^(1/2);
  Fa = rho*(mod_v^2)*Cd*A*(1/2)*(-x(7:9,k)/mod_v);

  x1_der=R*x(7:9,k);
  x2_der=Euler2Q(x(4:6,k))*x(10:12,k);
  x3_der=-skew(x(10:12,k))*x(7:9,k) + (1/m)*Fa + (1/m)*Fg +(1/m)*Fp(1:3,k);
  x4_der=inv(J)*(-skew(x(10:12,k)))*J*x(10:12,k) + inv(J)*np(1:3,k);


  fxu = [ x1_der;
          x2_der;#R *skew(x( 10:12 , k ));
          x3_der; %u(1:3) Força_Propulsão/Thrust + Força Gravitica, se indico o drone para impor um Thrust de 10N ([0 0 -10]N) o input será isso mais Fg(vetor) %u(2) é o momento np
          x4_der
        ];

  x(:,k+1) = x(:,k) + Dt*fxu ;
  y(:,k) = x(:,k) ; #Quero observar todos os estados
end

figure(90321);
plot(t,y(1:3,:));
grid on;
legend('$$p_x$$','$$p_y$$','$$p_z$$');

figure(90322);
plot(t,y(4:6,:));
grid on;
legend('$$fi$$','$$theta$$','$$psi$$');

figure(90323);
plot(t,y(7:9,:));
grid on;
legend('$$v_x$$','$$v_y$$','$$v_z$$');

figure(90324);
plot(t,y(10:12,:));
grid on;
legend('$$w_x$$','$$w_y$$','$$w_z$$');

figure(90325);
plot(t,u(1:4,:));
grid on;
legend('$$T_1$$','$$T_2$$','$$T_3$$','$$T_4$$');



