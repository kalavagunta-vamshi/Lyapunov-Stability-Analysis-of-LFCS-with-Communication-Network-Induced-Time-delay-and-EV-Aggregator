p1= 0.2;

ki= 7;

num= [0.1 0.7 1 0];

den= [1 7.08 10.56 20.8 7];

t= 0:.02:12;

c= -p1* step(num, den, t);
plot(t, c),xlabel('t, sec'), ylabel('pu') 
title('Frequency deviation step response')


