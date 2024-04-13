clear all 

J = 42.6e-6;
La = 170e-3;
Ra = 4.67;
B = 47.3e-6;
K = 14.7e-3;
Ka = 14.7e-3;

p3 = 1;
p2 = Ra/La + B/J;
p1 = (Ra*B+K*Ka)/(J*La);

tau = 0.5;
wc = 0:0.01:3.91;
for i = 1:size(wc,2)
    A1(:,i) = ((K)/(J*La))*wc(i)*sin(wc(i)*tau);
    B1(:,i) = ((K)/(J*La))*cos(wc(i)*tau);
    A2(:,i) = ((K)/(J*La))*wc(i)*cos(wc(i)*tau);
    B2(:,i) = -((K)/(J*La))*sin(wc(i)*tau);
    C1(:,i) = -p2*wc(i)^2;
    C2(:,i) = -p3*wc(i)^3+p1*wc(i);
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