clc
clear all
close all
%define values of Length, height, and width of cantilever beam
L=.20; %m
h=.009525; %m
b=.0254; %m
%defining all given constants in SI units
g=9.8; %m/s^2
E=3.4*10^9; %Pa
zeta=0.05;
H=0.00635; %m
D=.0127; %m
A_mag=3.981*10^-6;
b_mag=2.7;
rho=1180; %kg/m^3
m=.006; %kg
A_bar=b*h; %m^2
I=(b*h^3)/12; %m^4
%calculating a value used in transeq 1 and transeq2 from equation 12 of
%handout
msquig=m/(rho*A_bar*L);
%solve handout eq12 for p with beta=0
syms p
transeq1=1+cos(p*L)*cosh(p*L)+(msquig)*((p*L)*cos(p*L)*sinh(p*L)-p*L*sin(p*L)*cosh(p*L));
P_1=vpasolve(transeq1,p, [0 10]);
%calculate natural frequency of cantilever beam based on p value from
%transeq1
Omega_0=(P_1)^2*sqrt((E*I)/(rho*A_bar));
%calculate effective spring consant of cantilever beam
k_bar=3*E*I/(L^3);
%calculate critical distance from magnet to substrate before "pull in"
%happens
d_cr=(b_mag*A_mag/k_bar)^(1/(b_mag+1));
%calculate critical distance of shaking base, and deflection of beam at
%which "pull in" happens, choose and operating height slightly greater than
%the critical operating height
d_m_cr=H+d_cr+(1/k_bar)*(A_mag/(d_cr^(b_mag))+m*g);
delta_cr=d_m_cr-H-d_cr;
d_m_op=d_m_cr+.002;
%solve for distance from magnet to substrate, deflection of beam, and 
%effective spring constant due to magnetism at operating height
syms d_e_op
eq1=H-d_m_op+d_e_op+(1/k_bar)*(A_mag/((d_e_op)^b_mag));
D_E_OP=vpasolve(eq1, d_e_op);
delta_op=(1/k_bar)*(A_mag/D_E_OP^b_mag);
k_e_op=b_mag*A_mag/(D_E_OP^(b_mag+1));
%solve for p value of handout eq12 with non-zero Beta, now inlcudes
%magnetism effect
Beta=k_e_op*L^3/(E*I);
transeq2=1+cos(p*L)*cosh(p*L)+(msquig+Beta/((p*L)^4))*((p*L)*cos(p*L)*sinh(p*L)-p*L*sin(p*L)*cosh(p*L));
P_2=vpasolve(transeq2,p, [0 10]);
%calculate natural frequency from spring constants at operating height and
%a more accurate natural frequency based on solution to transeq 2 (handout
%eq 12 with magnetism included)
omega_m1=sqrt((k_bar-k_e_op)/m);
omega_m2=P_2^2*sqrt(E*I/(rho*A_bar));
%convert natural frequencies to Hertz from Radians per second
omega_m2_hz=omega_m2/(2*pi);
omega_0_hz=Omega_0/(2*pi);

%lines 67-94 repeat the calculations from lines 43 through 60 at 3 different operating
%heights, the 3 heights are the chosen height, and the chosen height plus
%or minus 0.5 mm. The amplitude of X0/Y0 eq4 in handout is then plotted
%against the ratio of chosen frequency to natural frequency for each
%operating height
syms d_e_op
matrix=[d_m_op d_m_op-.0005 d_m_op+.0005];

for i=1:3
    height=matrix(i);

    eq1=H-height+d_e_op+(1/k_bar)*(A_mag/((d_e_op)^b_mag));
    D_E_OP=vpasolve(eq1, d_e_op);
    delta_op=(1/k_bar)*(A_mag/D_E_OP^b_mag);

    k_e_op=b_mag*A_mag/(D_E_OP^(b_mag+1));
    Beta=k_e_op*L^3/(E*I);

    transeq2=1+cos(p*L)*cosh(p*L)+(msquig+Beta/((p*L)^4))*((p*L)*cos(p*L)*sinh(p*L)-p*L*sin(p*L)*cosh(p*L));
    P_2=vpasolve(transeq2,p, [0 10]);
    
    omega_m1=sqrt((k_bar-k_e_op)/m);
    omega_m2=P_2^2*sqrt(E*I/(rho*A_bar));
    
    omega_m2_hz=omega_m2/(2*pi);
    omega_0_hz=Omega_0/(2*pi);
    
    syms omega
    ratio=omega/omega_m2_hz;
    X0Y0m=(1+(2*zeta*ratio)^2)^(1/2)/((1-ratio^2)^2+(2*zeta*ratio)^2)^(1/2);
    fplot(X0Y0m, [0 102])
    hold on
end
title('Amplitude of Vibration at Different Frequencies')
xlabel('Operating Frequency (Hz)')
ylabel('X0/Y0')

fprintf('Chosen length is %s meters\n', L)
fprintf('Chosen base is %s meters\n', b)
fprintf('Chosen height is %s meters\n', h)
fprintf('Accurate non_magnetic natural frequency is %s radians per second\n', Omega_0)
fprintf('Effective spring constant of beam is %s Newtons per meter\n', k_bar)
fprintf('Critical gap, dc, is %s meters \n', d_cr)
fprintf('Critical beam, delta cr, deflection is %s meters\n', delta_cr)
fprintf('Critical base location, dmcr, is %s meters\n', d_m_cr)
fprintf('Chosen operating height, dm, is %s meters\n', d_m_op)
fprintf('Beam deflection at operating height, delta op, is %s meters\n', delta_op)
fprintf('Effective magnetic spring constant at operating height, ke op is %s Newtons per meter\n', k_e_op)
fprintf('Gap at operating height, de op, is %s meters\n', D_E_OP)
fprintf('Natural frequency at operating height with magnetic effect included is %s radians per second\n', omega_m1)
fprintf('More accurate natural frequency at operating height with magnetic effect included is %s radians per second\n', omega_m2)
