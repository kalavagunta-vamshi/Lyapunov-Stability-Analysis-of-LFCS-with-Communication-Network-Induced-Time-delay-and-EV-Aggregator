clear all
J = 42.6e-6;
L = 170e-3;
R = 4.67;
f = 47.3e-6;
Kt = 14.7e-3;
Kb = 14.7e-3;

Kp = 0.1;
Ki = 0.5;

p3 = 1;
p2 = R/L + f/J;
p1 = (R*f+Kt*Kb)/(J*L);
p0 = 0;

q1 = (Kt*Kp)/(J*L);
q0 = (Kt*Ki)/(J*L);

t6 = p3^2;
t4 = p2^2-2*p1*p3;
t2 = p1^2-q1^2;
t0 = -q0^2;
 
k=1;
r = roots([t6 t4 t2 t0]);
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
    w = wa(:,n);
    derW = 3*t6*w^4+2*t4*w^2+t2;
    RT = sign(derW);
    if(RT==1)
        m = m+1;
        P = p3*(j*w)^3+p2*(j*w)^2+p1*(j*w)^1;
        Q = q1*(j*w)^1+q0;
        SIN = imag(P/Q);
        COS = real(-P/Q);
        theta = atan(SIN/COS);
        if COS<0
            tau1(:,m) = ((theta+pi)/w);
        else
            tau1(:,m) = (theta/w);
        end
        RT;
        w;
        tau1;
    else
        RT;
        w;  
    end
end
Kp
Ki
format long g
tau=min(tau1)
format