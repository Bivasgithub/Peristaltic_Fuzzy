 function nanogeometry
l=0.1;
m=0.2;
a=0.1;
% z=input('Enter z: z=0:0.01:8;');  %z=2
z=0.2;
h=1+l*z+a*sin(2*pi*(z-m));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_1=400;  %%Cu
k_2=9.7;   %%Fe2O3 385
k_f=0.492;
alpha=0.15;
beta=0.15;
phi_1=0.1+0.05*alpha+beta*0.1-alpha*beta*0.1;
phi_2=0.1+0.05*alpha+beta*0.1-alpha*beta*0.1;
k_bf=k_f*((k_1+2*k_f-2*phi_1*(k_f-k_1))/(k_1+2*k_f+phi_1*(k_f-k_1)));
% k_hnf=k_f*(((k_1+2*k_f-2*phi_1*(k_f-k_1))/(k_1+2*k_f+phi_1*(k_f-...
%       k_1)))*((k_2*+2*k_bf-2*phi_2*(k_bf-k_2))/(k_2+2*k_bf+...
%       phi_2*(k_bf-k_2))));
A6=((k_1+2*k_f-2*phi_1*(k_f-k_1))/(k_1+2*k_f+phi_1*(k_f-...
   k_1)))*((k_2*+2*k_bf-2*phi_2*(k_bf-k_2))/(k_2+2*k_bf+phi_2*(k_bf-...
   k_2)));
% x1=1/(((1-phi_1)^2.5)*((1-phi_2)^2.5)*((1-phi_3)^2.5));
rho_1=8933;
rho_2=5180;  %10500
% rho_3=19320;
rho_f=1063;
bita_1=16.7;
bita_2=0.000013;    %18.7
% bita_3=14;
bita_f=1.8;
sigma_f=0.6670;
sigma_1=59600000;
sigma_2=25000;            %(6.3)*10^7
% sigma_3=(4.10)*10^7;
A4=(1-phi_2)*((1-phi_1)+phi_1*((rho_1*bita_1)/(rho_f*bita_f)))+...
   phi_2*((rho_2*bita_2)/(rho_f*bita_f));
sigma_bf=sigma_f*((sigma_1+2*sigma_f-2*phi_1*(sigma_f-...
         sigma_1))/(sigma_1+2*sigma_f+phi_1*(sigma_f-sigma_1)));
A3=sigma_bf*((sigma_2*+2*sigma_bf-2*phi_2*(sigma_bf-...
          sigma_2))/(sigma_2+2*sigma_bf+phi_2*(sigma_bf-sigma_2)));
A2=1/(((1-phi_1)^2.5)*((1-phi_2)^2.5));
%%%%%%%%%%%%%%%%%%%%%%%%temp%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W_s=0.1; %velocity slip;
% epsilon=0.01;  % radius of inner tube
% Q=0.1; %Heat sorce
% Rd=2.5;
% Br=1.5;
% Da=0.07;
% Fr=0.5;
% Re=0.01;
% M=0.5;%Hartmann number;
% Gr=0.1; %Grashoff number;
% lambda11=0.1;  
% eta=pi/2; 
W_s=0.05; %velocity slip;
epsilon=0.01;  % radius of inner tube
Q=0.1; %Heat sorce
Rd=0;
Br=0;
Da=0.09;
Fr=0.0;
Re=0.0;
M=2;%Hartmann number;
Gr=2; %Grashoff number;
lambda11=0.0;  
eta=pi/2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%velocity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W_s=0.05; %velocity slip;
% epsilon=0.1;  % radius of inner tube
% Q=1; %Heat sorce
% Rd=0.1;
% Br=1.5;
% Da=0.05;
% Fr=0.5;
% Re=0.95;
% M=1.5;%Hartmann number;
% Gr=1; %Grashoff number;
% lambda11=0.2;  
% eta=pi/2; 
%**************************************************************************
% % %TEMPERATURE PROFILE
% zpu=input('Enter z: z=0:0.01:8;');  %z=2
r=0:0.01:1;
for i=1:101
% h=NDStenosisGeometry(d,p,zpu);
a0=(A3*M^2*Br)/(A6+Rd);
a1=Q/(A6+Rd);
a2=(a0*(W_s)^2)/(h-epsilon)^2;
a3=(2*a0*((W_s)^2)*epsilon)/(h-epsilon)^2;
a4=(a0*((W_s)^2)*epsilon^2)/(h-epsilon)^2;
a5=(1+lambda11)/A2;
a6=((1+lambda11)/Da)+((A3*M^2*(1+lambda11))/A2);
a7=((Fr*Re*(1+lambda11))/A2);
a8=(A4*Gr*sin(eta)*(1+lambda11))/A2;
a9=W_s/(h-epsilon);
a10=a6*a9;
a11=a10*epsilon;
a12=a7*(a9^2);
a13=2*epsilon*a12;
a14=a12*(epsilon^2);
a15=a8/(epsilon-h);
a16=a15*h;
a17=-a11+a14+a16;
a18=a10-a13-a15;
D21=(1/log(epsilon/h))*(a9*(epsilon-h)-(a17/4)*(epsilon^2-h^2)-...
    (a18/9)*(epsilon^3-h^3)-(a12/16)*(epsilon^4-h^4));
D22=((epsilon^2-h^2)*(1+lambda11))/(4*A2*log(epsilon/h));
D11=-D21*log(h)+a9*h-(a17/4)*h^2-(a18/9)*h^3-(a12/16)*h^4;
D12=D22*log(h)-((h^2*(1+lambda11))/(4*A2));
a19=((a0*D11^2)/4)+((3*a0*D21^2)/8);
a191=((a0*D12^2)/4)+((3*a0*D22^2)/8);
C2=(1/log(epsilon/h))*(1+(a2/16)*(epsilon^4-h^4)-(a3/9)*(epsilon^3-h^3)+...
   ((a4+a1)/4)*(epsilon^2-h^2));
C1=-C2*log(h)-(h/(h-epsilon))+((a2*(h^4))/16)-((a3*(h^3))/9)+((a4+a1)/4)*(h^2);
C41=(1/(log(epsilon/h)))*(((a0*(D21^2))/4)*((epsilon*(log(epsilon)))^2-...
    (h*(log(h)))^2)-((a0*(D21^2))/2)*(((epsilon^2)*(log(epsilon)))-...
    ((h^2)*(log(h))))+a19*(epsilon^2-h^2)+...
    ((a0*(a9^2))/16)*(epsilon^4-h^4)+(a0/36)*((a17/4)^2)*(epsilon^6-h^6)+...
    (a0/64)*((a18/9)^2)*(epsilon^8-h^8)+(a0/100)*((a12/16)^2)*(epsilon^10-h^10));
C42=(1/(log(epsilon/h)))*((((epsilon^6-h^6)*a0*(1+lambda11)^2)/((24*A2)^2))+...
    (epsilon^2-h^2)*a191+(((a0*(D22^2))/4)*((epsilon*(log(epsilon)))^2-...
    (h*(log(h)))^2))-(((a0*(D22^2))/2)*(((epsilon^2)*(log(epsilon)))-...
    ((h^2)*(log(h))))));
C43=((a0*(D21^2))/4)*((h*(log(h)))^2)-((a0*(D21^2))/2)*(((h^2)*(log(h))))+...
    a19*h^2+((a0*(a9^2))/16)*(h^4)+(a0/36)*((a17/4)^2)*(h^6)+...
    (a0/64)*((a18/9)^2)*(h^8)+(a0/100)*((a12/16)^2)*(h^10);
C44=-C41*log(h)+C43;
C45=-C42*log(h)+a0*((1+lambda11)/(24*A2))^2*(h^6)+a191*(h^2)+...
     ((a0*(D22^2))/4)*(h*(log(h)))^2-((a0*(D22^2))/2)*((h^2)*(log(h)));
Lambda0=(h/(h-epsilon))+C1+C44;
Lambda1=C2+C41;
Lambda2=((a0*(D21^2))/4);
Lambda3=((a0*(D21^2))/2);
Lambda4=((a4+a1)/4)+a19;
Lambda5=(a2/16)+((a0*(a9^2))/16);
Lambda6=(a3/9);
Lambda7=(a0/36)*((a17/4)^2);
Lambda8=(a0/64)*((a18/9)^2);
Lambda9=(a0/100)*((a12/16)^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a20=a6*D11+a7*D11^2-a8*C1;
a21=a6*D21-a8*C2;
a22=a7*(D21^2);
a23=a6*a9-a15;
a24=((a6*a17)/4)+(a7*(a9^2))+((a8*(a4+a1))/4);
a25=(((a6*a18)-(a8*a3))/9);
a26=(((a6*a12)+(a7*(a17^2))+(a8*a2))/16);
a27=a7*(a18/9)^2;
a28=a7*(a12/16)^2;
a29=a6*D12;
a30=a6*D22;
a31=((a6*(1+lambda11))/(4*A2));
a32=a7*(D12^2);
a33=a7*(D22^2);
a34=((a7*((1+lambda11)^2))/((4*A2)^2));
a35=((a29+a30)/4);
a36=(((2*a32)+(3*a33))/8);
a37=(((2*a20)-(2*a21)+(3*a22))/8);
a38=((a21-(2*a22))/4);
D41=(1/(log(epsilon/h)))*(a37*(h^2-epsilon^2)+a38*(((h^2)*(log(h)))-...
    ((epsilon^2)*(log(epsilon))))+(a22/4)*((h*(log(h)))^2-...
    (epsilon*(log(epsilon)))^2)-(a23/9)*(h^3-epsilon^3)+...
    (a24/16)*(h^4-epsilon^4)+(a25/25)*(h^5-epsilon^5)+...
    (a26/36)*(h^6-epsilon^6)+(a27/64)*(h^8-epsilon^8)+(a28/100)*(h^10-...
    epsilon^10));
D42=-(1/(log(epsilon/h)))*(a35*(epsilon^2-h^2)-...
    (a30/4)*(((epsilon^2)*(log(epsilon)))-...
    ((h^2)*(log(h))))+(a31/16)*(epsilon^4-h^4));
D43=-(1/(log(epsilon/h)))*(a36*(epsilon^2-h^2)+(a34/36)*(epsilon^6-h^6)+...
    (a33/4)*((epsilon*(log(epsilon)))^2-(h*(log(h)))^2)-...
    (a33/2)*(((epsilon^2)*(log(epsilon)))-...
    ((h^2)*(log(h)))));
D31=-(a37*h^2+a38*h^2*(log(h))+(a22/4)*(h*(log(h)))^2-(a23/9)*h^3+...
     (a24/16)*h^4+(a25/25)*h^5+(a26/36)*h^6+(a27/64)*h^8+(a28/100)*h^10);
D32=a35*h^2-(a30/4)*h^2*(log(h))+(a31/16)*h^4;
D33=a36*h^2-(a33/2)*h^2*(log(h))+(a34/36)*h^6+(a33/4)*(h*(log(h)))^2;
D311=D31-D41*(log(h));
D312=D42*(log(h))+D32;
D313=D43*(log(h))+D33;
lambda1=((W_s*epsilon)/(epsilon-h))+D11+D311;
lambda2=D21+D41;
lambda3=a22/4;
lambda4=(a17/4)+a37;
lambda5=(a18-a23)/9;
lambda6=(a12+a24)/16;
lambda7=a25/25;
lambda8=a26/36;
lambda9=a27/64;
lambda10=a28/100;
lambda_11=D12-D312;
lambda12=D42-D22;
lambda13=a30/4;
lambda14=((1+lambda11)/(4*A2))+a35;
lambda15=a31/16;
B11=(lambda1/2)*(h^2-epsilon^2)+(lambda2/4)*(2*(((h^2)*(log(h)))-...
    ((epsilon^2)*(log(epsilon))))-(h^2-epsilon^2))+...
    (a38/16)*(4*(((h^4)*(log(h)))-((epsilon^4)*(log(epsilon))))-(h^4-...
    epsilon^4))+(lambda3/32)*(8*((((h^2)*(log(h)))^2)-...
    (((epsilon^2)*(log(epsilon)))^2))-4*(((h^4)*(log(h)))-...
    ((epsilon^4)*(log(epsilon))))+(h^4-epsilon^4))+(lambda4/4)*(h^4-...
    epsilon^4)+(lambda5/5)*(h^5-epsilon^5)+(lambda6/6)*(h^6-epsilon^6)+...
    (lambda7/7)*(h^7-epsilon^7)+(lambda8/8)*(h^8-epsilon^8)+...
    (lambda9/10)*(h^10-epsilon^10)+(lambda10/12)*(h^12-epsilon^12);
B12=(lambda_11/2)*(h^2-epsilon^2)+(lambda12/4)*(2*(((h^2)*(log(h)))-...
    ((epsilon^2)*(log(epsilon))))-(h^2-epsilon^2))-...
    (lambda13/16)*(4*(((h^4)*(log(h)))-((epsilon^4)*(log(epsilon))))-(h^4-...
     epsilon^4))+(lambda14/4)*(h^4-epsilon^4)+(lambda15/6)*(h^6-epsilon^6);
B13=(-D313/2)*(h^2-epsilon^2)+(D43/4)*(2*(((h^2)*(log(h)))-...
    ((epsilon^2)*(log(epsilon))))-(h^2-epsilon^2))+(a33/128)*(8*((((h^2)*(log(h)))^2)-...
    (((epsilon^2)*(log(epsilon)))^2))-4*(((h^4)*(log(h)))-...
    ((epsilon^4)*(log(epsilon))))+(h^4-epsilon^4))-(a33/2)*(4*(((h^4)*(log(h)))-...
    ((epsilon^4)*(log(epsilon))))-(h^4-...
     epsilon^4))+(a36/4)*(h^4-epsilon^4)+(a34/(36*8))*(h^8-epsilon^8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1 = 0.0;
E2 = 0.0;
E3 = 0.0;
% % pressure gradient define
dp_by_dz =  8 *pi^3 * a * ( - (E1 + E2) * cos(2 * pi * (z-m)) + (E3/(2 * pi)) * sin(2 * pi * (z-m)));
% F=0.1;
% dp_by_dz=((-B12+sqrt(B12^2-(4*B13*(B11-F))))/(2*B13));
% dp_by_dz=0.85;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Velocity profile%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w1(i)= lambda1+lambda2*log(r(i))+a38*(r(i))^2*log(r(i))+...
      lambda3*((r(i))*log(r(i)))^2+lambda4*(r(i))^2+lambda5*(r(i))^3+...
      lambda6*(r(i))^4+lambda7*(r(i))^5+lambda8*(r(i))^6+lambda9*(r(i))^8+...
      lambda10*(r(i))^10+(lambda_11+lambda12*log(r(i))-...
      lambda13*(r(i))^2*log(r(i))+lambda14*(r(i))^2+lambda15*(r(i))^4)*(dp_by_dz)+...
      (-D313+D43*log(r(i))+a36*(r(i))^2+(a34/36)*(r(i))^6+...
      (a33/4)*((r(i))*log(r(i)))^2-(a33/2)*(r(i))^2*log(r(i)))*(dp_by_dz)^2;
%%%%%% %%%%%%%%%%%%%%%%%%%Temperature profile%%%%%%%%%%%%%%%%%%%%%%%%%%
% T1(i)=Lambda0+Lambda1*log(r(i))-Lambda2*(r(i)*log(r(i)))^2+Lambda3*r(i)^2*log(r(i))-...
%    Lambda4*r(i)^2-Lambda5*r(i)^4+Lambda6*r(i)^3-Lambda7*r(i)^6-Lambda8*r(i)^8-Lambda9*r(i)^10+...
%    (C45+C42*log(r(i))-(((a0*(1+lambda11)^2)/((24*A2)^2))*r(i)^6+a191*r(i)^2+...
%    ((a0*(D22^2))/4)*(r(i)*(log(r(i))))^2-((a0*(D22^2))/2)*(r(i)^2*(log(r(i))))))*(dp_by_dz)^2;
end
% plot(r, w1, '-', 'LineWidth', 2);
% % axis([0.1 1 -2 0])
% % axis([0.1 0.88 -1 1])
% xlabel('r');
% ylabel('W');
% % ylabel('\tilde{\theta}');
% hold on;
% grid on;
Array1=csvread('Default Dataset (14).csv');
X1 = Array1(:, 1);
Y1 = Array1(:, 2);
plot(X1, Y1,'--g','LineWidth',2)
xlabel('r');
ylabel('u^*');
hold on;
Array1=csvread('Dataset (9).csv');
X1 = Array1(:, 1);
Y1 = Array1(:, 2);
plot(X1, Y1,'--b','LineWidth',2)
xlabel('r');
ylabel('u^*');
hold on;
end
