tic

clc
clf
close all
clear
% _____________________________________________

% E_electron=3000;  % MeV
% Mass_electron=0.5109989461;  % MeV
% gamma=E_electron/Mass_electron;   % Lorentz factor

gamma=7460.52;
beta=sqrt(1-gamma^(-2));    % Beta factor
v=beta*299792458;     % Beta*Light speed   (m/s)   &&  Light speed=299792458 m/s

f=logspace(2,14,100);     % Frequency  (Hz)
omega=2*pi*f;
t=5E-6/v:5E-6/v:1E-2/v;

C= 528;  % Circumference    (m)
T=C/v;   % Revolution period (s)   1/f   
f_rev=1/T;   % Revolution frequency  (Hz)  C/v
omega_rev=2*pi*f_rev;

n=omega/omega_rev;

% _____________________________________________
sigma1=1/(9.1E-07);      % NEG conductivity (S/m)
sigma2=1/(1.68E-08);  % copper conductivity (S/m)

tau1=0;   
tau2=0;         % Copper AC relaxation time (s)
% _____________________________________________
delta=2E-6;

b2=10E-3+delta; 
b3=11E-3+delta;  % outer dimeter (m)
% _____________________________________________


[Z1,Z2,Z3]=Impedance(f,gamma,sigma1,sigma2,tau1,tau2,delta,b2,b3);

W1=DiscreteFastFourierTransform(Z1,t,omega);
W2=DiscreteFastFourierTransform(Z2,t,omega);
W3=DiscreteFastFourierTransform(Z3,t,omega);

figure
semilogx(f,real(Z1),'-.')
hold on
semilogx(f,imag(Z1),'--')
grid on
title('$ Longitudinal Impedance $','Interpreter','latex')
ylabel('$ Z_{\parallel} \left( \Omega \right) $','Interpreter','latex')
xlabel('$ \omega \left( Hz \right) $','Interpreter','latex')

legend('$\Re{(Z_{\parallel})}$','$\Im{(Z_{\parallel})}$','Interpreter','latex','Location','northwest')

% load ZlongWSLSII.txt
% semilogx(ZlongWSLSII(:,1),ZlongWSLSII(:,2))
% semilogx(ZlongWSLSII(:,1),ZlongWSLSII(:,3))
figure
semilogx(f,real(Z1./n))
hold on

semilogx(f,imag(Z1./n))

grid on


figure
semilogx(f,real(Z2),'-.')
hold on
semilogx(f,imag(Z2),'--')
grid on
title('$ Dipolar Transverse Impedance $','Interpreter','latex')
ylabel('$ Z_{x}^{Dipolar} \left( \Omega /m \right) $','Interpreter','latex')
xlabel('$ \omega \left( Hz \right) $','Interpreter','latex')
legend('$\Re{(Z_{x}^{Dipolar})}$','$\Im{(Z_{x}^{Dipolar})}$','Interpreter','latex')

% load ZxdipWSLSII.txt
% semilogx(ZxdipWSLSII(:,1),ZxdipWSLSII(:,2))
% semilogx(ZxdipWSLSII(:,1),ZxdipWSLSII(:,3))

figure
semilogx(f,real(Z3),'-.')
hold on
semilogx(f,imag(Z3),'--')
grid on
title('$ Quadrupolar Transverse Impedance $','Interpreter','latex')
ylabel('$ Z_{x}^{Quadrupole} \left( \Omega /m \right) $','Interpreter','latex')
xlabel('$ \omega \left( Hz \right) $','Interpreter','latex')
legend('$\Re{(Z_{x}^{Quadrupole})}$','$\Im{(Z_{x-Quadrupole})}$','Interpreter','latex','Location','northwest')



% load ZxquadWSLSII.txt
% semilogx(ZxquadWSLSII(:,1),ZxquadWSLSII(:,2))
% semilogx(ZxquadWSLSII(:,1),ZxquadWSLSII(:,3))

figure
semilogx(v*t,real(W1)/pi)
load WlongWSLSII.txt
hold on
semilogx(WlongWSLSII(:,1),WlongWSLSII(:,2))

grid on

figure
semilogx(v*t,imag(W2)/pi)
load WxdipWSLSII.txt
hold on
semilogx(WxdipWSLSII(:,1),WxdipWSLSII(:,2))
grid on


figure
semilogx(v*t,imag(W3)/pi)
load WxquadWSLSII.txt
hold on
semilogx(WxquadWSLSII(:,1),WxquadWSLSII(:,2))
grid on



toc
