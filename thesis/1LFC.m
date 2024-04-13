Pl= 0.2; num= [0.1 0.7 1];
den= [1 7.08 10.56 20.8];

t= 0:.02:10;

c= -pl* step(num,den,t);

plot(t, c), xlabel(‘t, sec’), ylabel(‘pu’) title(‘Frequency deviation step response’),grid timespec(num, den)
