p0=4443.577691;
Ta0=273;
Ra=287;
m1=0.0065;
H1=11000;
H2=20000;
m2=0.002;
g=9.81

%for h >H2
CC=Ta0-m1*H1-m2*H2;
A=Ra*CC/g;
B=Ra*m2/g;
[H,P] = ode15s(@myfunction3,[20000,100000],p0);

%for h > H1 && h<H2
p0=20453.68119;
C=g/(Ra*(Ta0-m1*H1));
[H,P] = ode15s(@myfunction2,[11000,20000],p0);

%for h <H1
p0=101000;
D=-Ra*m1/g;
E=Ta0*Ra/g;
[H,P]=ode15s(@myfunction1,[0,11000],p0);