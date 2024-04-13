clear all

M=8.8;
D=1;
Tg=0.2;
Tc=0.3;
Tr=12;
Fp=1/6;
R=1/11;
beta=21;
Kev=1;
Tev=0.1;
 
alpha1=0;
alpha0=1;
 
Q4=[alpha1*beta*R*Kev*Tg*Tr*Tc];
Q3=[alpha1*beta*R*Kev*(Tr*Tc+Tg*Tc+Tg*Tr)];
Qq3=[alpha1*beta*R*Kev*Tg*Tr*Tc];
Q2=[alpha1*beta*R*Kev*(Tc+Tr+Tg)];
Qq2=[alpha1*beta*R*Kev*(Tr*Tc+Tc*Tg*Tc+Tg*Tr)];
Qq0=[alpha1*beta*R*Kev];
Q1=[alpha1*beta*R*Kev];
Qq1=[alpha1*beta*R*Kev*(Tc+Tr+Tg)];
 
p6=M*R*Tg*Tr*Tc*Tev;
p5=D*R*Tg*Tc*Tr*Tev+M*R*[Tg*Tr*Tc+Tr*Tc*Tev+Tg*Tc*Tev+Tg*Tr*Tev];
p4=D*R*[Tr*Tg*Tc+Tr*Tc*Tev+Tg*Tc*Tev+Tg*Tr*Tev]+M*R*[Tr*Tc+Tg*Tc+Tg*Tr+Tc*Tev+Tr*Tev+Tg*Tev];
p3=D*R*[Tr*Tc+Tg*Tc+Tg*Tr+Tc*Tev+Tr*Tev+Tg*Tev]+M*R*[Tc+Tr+Tg+Tev]+Fp*Tr*Tev;
P3=[alpha0*beta*R*Fp*Tr*Tev];
p2=D*R*[Tc+Tr+Tg+Tev]+M*R+Fp*Tr+Tev;
P2=[alpha0*beta*R*[Tev+Fp*Tr]];
Pp2=[alpha0*beta*R*Fp*Tr*Tev];
p1=D*R+1;
P1=[alpha0*beta*R+1];
Pp1=[alpha0*beta*R*[Tev+Fp*Tr]];
Pp0=[alpha0*beta*R];
 
t = 0.5;
wc = 0:0.001:8;
for i=1:size(wc,2)
    A1(:,i) = [-P2*wc(i)^2]+[Q4*wc(i)^2*cos(wc(i)*t)]-[Q2*wc(i)^2*cos(wc(i)*t)]-[Q3*wc(i)^3*sin(wc(i)*t)]+[Q1*wc(i)*sin(wc(i)*t)];
    B1(:,i) = [-Pp2*wc(i)^2]+[Pp0-Qq2*wc(i)^2*cos(wc(i)*t)]+[Qq0*cos(wc(i)*t)]-[Qq3*wc(i)^3*sin(wc(i)*t)]+[Qq1*wc(i)*sin(wc(i)*t)];
    A2(:,i) = [-P3*wc(i)^3]+[P1*wc(i)]-[Q3*wc(i)^3*cos(wc(i)*t)]+[Q4*wc(i)^4*sin(wc(i)*t)]+[Q2*wc(i)^2*sin(wc(i)*t)];
    B2(:,i) = [Pp1*wc(i)]-[Qq3*wc(i)^3*cos(wc(i)*t)]+[Qq1*wc(i)*cos(wc(i)*t)]+[Qq2*wc(i)^2*sin(wc(i)*t)]-[Qq0*sin(wc(i)*t)];
    C1(:,i) = [-p6*wc(i)^6]+[p4*wc(i)^4]-[p2*wc(i)^2];
    C2(:,i) = [p5*wc(i)^5]-[p3*wc(i)^3]+[p1*wc(i)];
    Kp(:,i) = (B1(:,i)*C2(:,i)-B2(:,i)*C1(:,i))/(A1(:,i)*B2(:,i)-A2(:,i)*B1(:,i));
    Ki(:,i) = (A2(:,i)*C1(:,i)-A1(:,i)*C2(:,i))/(A1(:,i)*B2(:,i)-A2(:,i)*B1(:,i));
end
plot(Kp,Ki)
hold on
 
 
y =-1:0.001:18;
x = 0.001*y;
plot(y,x)
hold on
grid
axis([-2 18 -1 8])
xlabel('Proportional Controller Gain K_{P}')
ylabel('Integral Controller Gain K_{I}')
title('Stability Region for \tau=0.5')
 
 
 
 
 
 


