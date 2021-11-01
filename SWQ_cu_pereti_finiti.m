% Valori si functii proprii ale hamiltonianului in groapa cuantica
% dreptunghiulara (SQW) cu pereti finiti
clc; clear; close all;
h=6.626*1e-34; % J*s; constanta Planck
hbar=h/2/pi; % constanta Planck redusa
eV=1.602*1e-19; % J; electron-Volt (unitate atomica de energie)
m0=9.1*1e-31; % kg; masa electronului in vid
meff=0.067*m0; % masa efectiva a electronului de conductie in GaAs
U0=0.280*eV; % bariera de energie potentiala GaAs/AlGaAs
w=12*1e-9; % m; largimea gropii de potential
xs=-1.5*w; xd=1.5*w; Nx=1000; x=linspace(xs,xd,Nx); dx=x(2)-x(1); % discretizarea x
a=2*meff*dx^2/hbar^2; % vezi Curs 5 partea a II-a, rel. (***)
U=zeros(1,Nx); U(x<-w/2)=U0; U(x>w/2)=U0; % profilul de energie potentiala
phi=zeros(1,Nx); % prealocare functie de unda
phi(2)=1e-6; % pentru a evita o solutie identic nula
Emin=0; Emax=U0; NE=5000; E=U0*logspace(-2,0,NE); % progresie exponentiala
phid=zeros(1,NE);
for iE=1:NE % stabileste energia de tir
    for ix=2:Nx-1 % implementeaza recurenta Schrodinger
        phi(ix+1)=2*phi(ix)-phi(ix-1)+a*(U(ix)-E(iE))*phi(ix); % Curs 5
    end
    phid(iE)=phi(Nx); % valoarea functiei de unda la frontiera dreapta
end
logphi=log(abs(phid)); % logaritmare pentru postselectie mai rapida
contor=1; Ep=zeros(1,1);
for iE=2:NE-1
    if (logphi(iE-1)>logphi(iE))&&(logphi(iE)<logphi(iE+1)) % test minim
        Ep(contor)=E(iE); contor=contor+1;
    end
end
figure(1);
plot(E/eV*1000,logphi,'-r');
xlabel('E tir / meV'); ylabel('ln(abs(phid))'); grid;
title('Curba de selectie a energiilor proprii');
disp('Valorile proprii ale energiei (meV):'); disp(Ep/eV*1000)
figure(2);
NEp=length(Ep); % numarul starilor proprii
plot(x/1e-9,U/eV*1000,'-r'); hold on % profilul energiei potentiale
xlabel('x/nm'); ylabel('U/meV, phi'); grid;
title('Nivele de energie si functii proprii');
for iEp=1:NEp % ciclul energiilor proprii
   for ix=2:Nx-1 % implementeaza recurenta Schrodinger
     phi(ix+1)=2*phi(ix)-phi(ix-1)+a*(U(ix)-Ep(iEp))*phi(ix);
   end
   plot(x/1e-9,Ep(iEp)/eV*1000*ones(1,Nx),'-g'); % nivel de energie
   plot(x/1e-9,Ep(iEp)/eV*1000+phi/max(phi)*U0/10/eV*1000,'-b');
end
axis([-1.25*w/1e-9,1.25*w/1e-9,-0.1*U0/eV*1000,1.1*U0/eV*1000])