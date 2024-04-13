Clear all
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
alpha1=0.2;
alpha0=0.8;
Q4=Kp*[alpha1*beta*R*Kev*Tg*Tr*Tc];
Q3= Kp*[alpha1*beta*R*Kev*[Tr*Tc+Tg*Tc+Tg*Tr]];
Qq3=Ki*[alpha1*beta*R*Kev*Tg*Tr*Tc];
Q2=Kp*[alpha1*beta*R*Kev*[Tc+Tr+Tg]];
Qq2=Ki*[alpha1*beta*R*Kev*[Tr*Tc+Tc*Tg*Tc+Tg*Tr];
Qq0=Ki*[alpha1*beta*R*Kev*Ki];
Q1=Kp*[alpha1*beta*R*Kev];
Qq1=Ki*[alpha1*beta*R*Kev*[Tc+Tr+Tg]];
p6=M*R*Tg*Tr*Tc*Tev;
p5=D*R*Tg*Tc*Tr*Tev+M*R*[Tg*Tr*Tc+Tr*Tc*Tev+Tg*Tc*Tev+Tg*Tr*Tev];
p4=D*R*[Tr*Tg*Tc+Tr*Tc*Tev+Tg*Tc*Tev+Tg*Tr*Tev]+M*R*[Tr*Tc+Tg*Tc+Tg*Tr+Tc*Tev+Tr*Tev+Tg*Tev];
p3=D*R*[Tr*Tc+Tg*Tc+Tg*Tr+Tc*Tev+Tr*Tev+Tg*Tev]+M*R*[Tc+Tr+Tg+Tev]+Fp*Tr*Tev;
P3=Kp*[alpha0*beta*R*Fp*Tr*Tev];
p2=D*R*[Tc+Tr+Tg+Tev]+M*R+Fp*Tr+Tev;
P2=Ki*[alpha0*beta*R*[Tev+Fp*Tr]];
Pp2=Ki*[alpha0*beta*R*Fp*Tr*Tev];
p1=D*R+1;
P1=Kp*[alpha0*beta*R+1];
Pp1=Ki*[alpha0*beta*R*[Tev+Fp*Tr]];
Pp0=Ki*[alpha0*beta*R];
t = 0.5;
wc = 0:0.01:3.91;
for i=1:size(wc,2)
    A1(:,i) = -P2*wc^2+Q4*wc^2*cos(wc*t)-Q2*wc^2*cos(wc*t)-Q3*wc^3*sin(wc*t)+Q1*wc*sin(wc*t);
    B1(:,i) = -Pp2*wc^2+Pp0-Qq2*wc^2*cos(wc*t)+Qq0*cos(wc*t)-Qq3*wc^3*sin(wc*t)+Qq1*wc*sin(wc*t);
    A2(:,i) = -P3*wc^3+P1*wc-Q3*wc^3*cos(wc*t)+Q4*wc^4*sin(wc*t)+Q2*wc^2*sin(wc*t);
    B2(:,i) = Pp1*wc-Qq3*wc^3*cos(wc*t)+Qq1*wc*cos(wc*t)+Qq2*wc^2*sin(wc*t)-Qq0*sin(wc*t);
    C1(:,i) = -p6*wc^6+p4*wc^4-p2*wc^2;
    C2(:,i) = p5*wc^5-p3*wc^3+p1*wc;
    Kp(:,i) = (B1(:,i)*C2(:,i)-B2(:,i)*C1(:,i))/(A1(:,i)*B2(:,i)-A2(:,i)*B1(:,i));
    Ki(:,i) = (A2(:,i)*C1(:,i)-A1(:,i)*C2(:,i))/(A1(:,i)*B2(:,i)-A2(:,i)*B1(:,i));
end
plot(Kp,Ki)
hold on

x = 0;
y = -0.03:0.0001:0.059;
plot(y,x,'k')
grid
axis([-0.04 0.07 -0.02 0.12])
xlabel('Proportional Controller Gain K_{P}')
ylabel('Integral Controller Gain K_{I}')
title('Stability Region for \tau=0.5')

