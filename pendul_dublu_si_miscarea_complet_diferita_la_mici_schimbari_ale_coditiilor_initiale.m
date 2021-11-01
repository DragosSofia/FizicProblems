% Pendulul gravitational dublu

clear; close all; clc;
RT=1; % selector real time(1) / slow motion(0)
g=9.80665; % m/s^2; acceleratia gravitationala terestra
% Parametrii fizici ai sistemului mecanic:
L1=2.2; L2=1.3; % m; lungimile tijelor
m1=1.2; m2=1.7; % kg; masele corpurilor
% Conditiile initiale - orice unghiuri intre -180 si +180!
theta10=100; theta20=-10; % grade; unghiurile initiale
theta10=theta10*pi/180; theta20=theta20*pi/180; % conversia in rad
OM10=-30; OM20=20; % grade/s; vitezele unghiulare initiale
OM10=OM10*pi/180; OM20=OM20*pi/180; % conversia in rad/s
% Definirea duratelor caracteristice:
omega1=sqrt(g/L1); omega2=sqrt(g/L2); % pulsatii proprii ale componentelor
T1=2*pi/omega1; T2=2*pi/omega2; % perioade proprii ale componentelor
T=max(T1,T2); % "timp caracteristic" al miscarii pendulului dublu
ti=0; tf=10*T; N=20000; t=linspace(ti,tf,N); dt=t(2)-t(1); % timpul discret
% Prealocare si valori de start petru primul corp:
theta11=zeros(1,N); theta21=theta11; % prealocare unghiuri
OM11=zeros(1,N); OM21=OM11; % prealocare viteze unghiulare
theta11(1)=theta10; theta21(1)=theta20; % unghiuri de start pas 1 al corpului 1
theta11(2)=theta10+OM10*dt; theta21(2)=theta20+OM20*dt; % unghiuri de start pas 2
OM11(1)=OM10; OM21(1)=OM20; % valori de start ale vitezelor unghiulare
% Notatii ajutatoare:
miu=1+m1/m2; % coeficient adimensional
r=L2/L1; % coeficient adimensional
a11=miu; a22=r; % coeficienti diagonala principala (constanti)
tic;
for i=2:N-1 % ciclul recurentelor
    aux=theta21(i)-theta11(i);
    a21=cos(aux); a12=a21*r; % coeficienti diagonala secundara (variabili)
    % Pentru vitezele unghiulare curente folosim derivatele la stanga:
    OM12(i)=(theta11(i)-theta11(i-1))/dt; % viteza corpului 1 la pasul i
    OM22(i)=(theta21(i)-theta21(i-1))/dt; % viteza corpului 2 la pasul i
    b1=r*OM21(i)^2*sin(aux)-g/L1*miu*sin(theta11(i)); % termen "liber" 1
    b2=-OM11(i)^2*sin(aux)-g/L1*sin(theta21(i)); % termen "liber" 2
    A=[a11,a12;a21,a22]; B=[b1;b2]; % matrice sistem si coloana termeni liberi
    E=A\B; % rezolvarea sistemului liniar in forma matriceala
    eps1=E(1); eps2=E(2); % acceleratiile unghiulare curente
    % Recurentele de ordinul II:
    theta11(i+1)=2*theta11(i)-theta11(i-1)+dt^2*eps1; % corp 1
    theta21(i+1)=2*theta21(i)-theta21(i-1)+dt^2*eps2; % corp 2
end;

toc; % afiseaza timpul de calcul al solutiei numerice
% Coordonate carteziene ale corpurilor:
x11=L1*sin(theta11); x21=x11+L2*sin(theta21); % coordonate orizontale
y11=-L1*cos(theta11); y21=y11-L2*cos(theta21); % coordonate verticale

% Prealocare si valori de start petru cel de al doilea corp:
theta12=zeros(1,N); theta22=theta12; % prealocare unghiuri
OM12=zeros(1,N); OM22=OM12; % prealocare viteze unghiulare
theta12(1)=theta10; theta22(1)=theta20; % unghiuri de start pas 1 al corpului 1
theta12(2)=theta10 + ( OM10 + 0.001 )*dt; theta22(2)=theta20+( OM20 + 0.001 )*dt; % unghiuri de start pas 2
OM12(1)= OM10 + 0.001; OM22(1)= OM20 + 0.001; % valori de start ale vitezelor unghiulare

for i=2:N-1 % ciclul recurentelor
    aux=theta22(i)-theta12(i);
    a21=cos(aux); a12=a21*r; % coeficienti diagonala secundara (variabili)
    % Pentru vitezele unghiulare curente folosim derivatele la stanga:
    OM12(i)=(theta12(i)-theta12(i-1))/dt; % viteza corpului 1 la pasul i
    OM22(i)=(theta22(i)-theta22(i-1))/dt; % viteza corpului 2 la pasul i
    b1=r*OM22(i)^2*sin(aux)-g/L1*miu*sin(theta12(i)); % termen "liber" 1
    b2=-OM12(i)^2*sin(aux)-g/L1*sin(theta22(i)); % termen "liber" 2
    A=[a11,a12;a21,a22]; B=[b1;b2]; % matrice sistem si coloana termeni liberi
    E=A\B; % rezolvarea sistemului liniar in forma matriceala
    eps1=E(1); eps2=E(2); % acceleratiile unghiulare curente
    % Recurentele de ordinul II:
    theta12(i+1)=2*theta12(i)-theta12(i-1)+dt^2*eps1; % corp 1
    theta22(i+1)=2*theta22(i)-theta22(i-1)+dt^2*eps2; % corp 2
end;

toc; % afiseaza timpul de calcul al solutiei numerice
% Coordonate carteziene ale corpurilor:
x12=L1*sin(theta12); x22=x12+L2*sin(theta22); % coordonate orizontale
y12=-L1*cos(theta12); y22=y12-L2*cos(theta22); % coordonate verticale

% Prealocare si valori de start petru cel de al treilea corp:
theta13=zeros(1,N); theta23=theta13; % prealocare unghiuri
OM13=zeros(1,N); OM23=OM13; % prealocare viteze unghiulare
theta13(1)=theta10; theta23(1)=theta20; % unghiuri de start pas 1 al corpului 1
theta13(2)=theta10 + ( OM10 + 0.0001 )*dt; theta23(2)=theta20+( OM20 + 0.0001 )*dt; % unghiuri de start pas 2
OM13(1)= OM10 + 0.0001; OM22(1)= OM20 + 0.0001; % valori de start ale vitezelor unghiulare

for i=2:N-1 % ciclul recurentelor
    aux=theta23(i)-theta13(i);
    a21=cos(aux); a12=a21*r; % coeficienti diagonala secundara (variabili)
    % Pentru vitezele unghiulare curente folosim derivatele la stanga:
    OM13(i)=(theta13(i)-theta13(i-1))/dt; % viteza corpului 1 la pasul i
    OM23(i)=(theta23(i)-theta23(i-1))/dt; % viteza corpului 2 la pasul i
    b1=r*OM23(i)^2*sin(aux)-g/L1*miu*sin(theta13(i)); % termen "liber" 1
    b2=-OM13(i)^2*sin(aux)-g/L1*sin(theta23(i)); % termen "liber" 2
    A=[a11,a12;a21,a22]; B=[b1;b2]; % matrice sistem si coloana termeni liberi
    E=A\B; % rezolvarea sistemului liniar in forma matriceala
    eps1=E(1); eps2=E(2); % acceleratiile unghiulare curente
    % Recurentele de ordinul II:
    theta13(i+1)=2*theta13(i)-theta13(i-1)+dt^2*eps1; % corp 1
    theta23(i+1)=2*theta23(i)-theta23(i-1)+dt^2*eps2; % corp 2
end;

toc; % afiseaza timpul de calcul al solutiei numerice
% Coordonate carteziene ale corpurilor:
x13=L1*sin(theta13); x23=x13+L2*sin(theta23); % coordonate orizontale
y13=-L1*cos(theta13); y23=y13-L2*cos(theta23); % coordonate verticale

figure(1);
Lmax=L1+L2; % semilatura cadrului grafic
coef=30; % controleaza dimensiunile grafice ale corpurilor
rg1=coef*m1^(1/3); rg2=coef*m2^(1/3); % raze "grafice"
tic; simt=0; % porneste cronometrul si initializeaza timpul simularii
while simt<=tf % ciclul grafic
  j=abs(t-simt)==min(abs(t-simt)); % cauta cel mai apropiat t din discretizare
  %afisarea primului corp
  plot([0 x11(j) x21(j)],[0 y11(j) y21(j)],'-g','LineWidth',3); hold on; % tije
  xlabel('x/m'); ylabel('y/m');
  plot(0,0,'.k','MarkerSize',10); % articulatie de suspensie
  plot(x11(j),y11(j),'.r','MarkerSize',rg1); % corpul 1
  plot(x11(j),y11(j),'.k','MarkerSize',10); % articulatie corp 1
  plot(x21(j),y21(j),'.b','MarkerSize',rg2); % corpul 2
 
  
  %afisarea celui de al doilea corp
  plot([0 x12(j) x22(j)],[0 y12(j) y22(j)],'-g','LineWidth',3); hold on; % tije
  xlabel('x/m'); ylabel('y/m');
  plot(0,0,'.k','MarkerSize',10); % articulatie de suspensie
  plot(x12(j),y12(j),'.r','MarkerSize',rg1); % corpul 1
  plot(x12(j),y12(j),'.k','MarkerSize',10); % articulatie corp 1
  plot(x22(j),y22(j),'.b','MarkerSize',rg2); % corpul 2
  
  %afisarea celui de al treilea corp
  plot([0 x13(j) x23(j)],[0 y13(j) y23(j)],'-g','LineWidth',3); hold on; % tije
  xlabel('x/m'); ylabel('y/m');
  plot(0,0,'.k','MarkerSize',10); % articulatie de suspensie
  plot(x13(j),y13(j),'.r','MarkerSize',rg1); % corpul 1
  plot(x13(j),y13(j),'.k','MarkerSize',10); % articulatie corp 1
  plot(x23(j),y23(j),'.b','MarkerSize',rg2); % corpul 2
  
  axis([-Lmax Lmax -Lmax Lmax]); axis square; % cadrul grafic
  if RT==1 % real time(1) / slow motion(0)
    simt=toc; % actualizeaza timpul simularii cu ceasul sistemului
    text(3/5*Lmax,4/5*Lmax,['t = ',num2str(round(t(j))),' s']);
  else
    simt=simt+1e-2; % incrementeaza cu o centisecunda
    text(3/5*Lmax,4/5*Lmax,['t=',num2str(round(t(j)*100)),' cs']);
  end
  pause(1e-6); hold off
end