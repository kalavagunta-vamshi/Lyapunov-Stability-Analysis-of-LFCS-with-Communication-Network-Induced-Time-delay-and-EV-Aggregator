clear ALL
 
M=    8.8;
D=    1;
Tg=   0.2;
Tc=   0.3;
Tr=   12;
Fp=   1/6;
R=    1/11;
Beta= 21;
Kev=  1;
Tev=  0.1;
 
alpha0= 0.8;
alpha1= 0.2;
 
Ki = 0.6; 
Kp = 0.2;
 
p6= (M*R*Tg*Tr*Tc*Tev);
p5= (D*R*Tg*Tr*Tc*Tev) + (M*R)*(Tg*Tr*Tc + Tr*Tc*Tev + Tg*Tc*Tev +  Tg*Tr*Tev);
p4= (D*R)*(Tg*Tr*Tc + Tr*Tc*Tev + Tg*Tc*Tev + Tg*Tr*Tev) + (M*R)*(Tr*Tc + Tg*Tc + Tg*Tr + Tc*Tev + Tg*Tev+ Tr*Tev); 
p3= (D*R)*( Tc*Tr+ Tg*Tc + Tg*Tc + Tg*Tr + Tc*Tev + Tr*Tev + Tg*Tev) + (M*R)*(Tc+Tr+Tg+Tev)+(Fp*Tr*Tev) + (alpha0*Beta*R*Kp*Fp*Tr*Tev);
p2= (D*R)*(Tc+Tr+Tg+Tev)+(M*R)+(Fp*Tr)+Tev+(alpha0*Beta*R)*(Kp*Tev + Kp*Fp*Tr + Ki*Fp*Tr*Tev);
p1= (D*R) + 1 + (alpha0*Beta*R)*(Kp+ Ki*Tev + Ki*Fp*Tr);
p0= (alpha0*Beta*R*Ki);
 
 
q4= (alpha1*Beta*R*Kev*Kp*Tg*Tr*Tc);
q3= (alpha1*Beta*R*Kev)*(Kp*Tr*Tc + Kp*Tg*Tc + Kp*Tg*Tr + Ki*Tg*Tc*Tr);
q2= (alpha1*Beta*R*Kev)*(Kp*Tc + Kp*Tr + Kp*Tg + Ki*Tr*Tc + Ki*Tg*Tc +Ki*Tg*Tr);
q1= (alpha1*Beta*R*Kev)*(Kp + Ki*Tc + Ki*Tr + Ki*Tg);
q0= alpha1*Beta*R*Kev*Ki;
 
 
t12=  p6^2;
t10=  p5^2 - 2*p6*p4;
t8 =  p4^2 + 2*p6*p2 - 2*p5*p3 - q4^2;
t6 =  p3^2 - 2*p6*p0 - 2*p4*p2 + 2*p5*p1 + 2*q4*q2 - q3^2;
t4 =  p2^2 + 2*p4*p0 - 2*p3*p1 - 2*q4*q0 + 2*q3*q1 - q2^2;
t2 =  p1^2 - 2*p2*p0 + 2*q2*q0 - q1^2;
t0 =  p0^2 - q0^2;
 
k=1;
r = roots([t12 t10 t8 t6 t4 t2 t0]);
for i=1:size(r,1)
    an(i)=angle(r(i))*(180/pi);
    if (an(i)==0)
        w1(k)=r(i);
        k=k+1;
    end
end
 
m = 0;
for n = 1:size(w1,2)
    wa(:,n) = sqrt(w1(:,n));
    wc = wa(:,n);
    derWc = 6*t12*wc^10 + 5*t10*wc^8 + 4*t8*wc^6 + 3*t6*wc^4 + 2*t4*wc^2 + t2;
    RT = sign(derWc);
    if(RT==1)
        m = m+1;
        P = p6*(j*wc)^6 + p5*(j*wc)^5 + p4*(j*wc)^4 + p3*(j*wc)^3 + p2*(j*wc)^2 + p1*(j*wc)^1 + p0;
        Q = q4*(j*wc)^4 + q3*(j*wc)^3 + q2*(j*wc)^2 + q1*(j*wc)^1 + q0;
        SIN = imag(P/Q);
        COS = real(-P/Q);
        theta = atan(SIN/COS);
        if COS<0
            tau1(:,m) = ((theta+pi)/wc);
        else
            tau1(:,m) = (theta/wc);
        end
        RT;
        wc;
        tau1;
    else
        RT;
        wc;  
    end
end
Kp
Ki
format long g
tau =  min(tau1)
format
 
 


