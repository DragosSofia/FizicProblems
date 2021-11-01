clc; clear; close all; %linie de igena

Gv = 1;%componenta vitezei in functie de timp
Gc = 1;% legile de miscare
T = 1; %traietoria 
F = 1; % reprezentarea dinamica

% parametri fizici
g = 9.80665; % acceleratia gravitationala standard [N/kg]
m = 1.2; % masa in proiectilului [kg]
G = m*g ; % forta de greutate [N]

% conditii initiale
v0 = 500; % viteza initiala [m/s] 
alpha0 = 43 ; % unghiul de lansare [grade]

% valori plauzibile ale coeficientilor de frecare fluida:
r1 = G / 2 / v0 ; % coeficientul teremnului liniar ar f de frecare
r2 = G / 2 / v0^2 ; % coeficientul teremnului patratic ar f de frecare

%Definirea intervalului de timp de interes 
t0= 0 ; tf = 2*v0/g*sind( alpha0 ); % timpul initial si supraestimarea timpului de zbor
N = 1000 ; t = linspace( t0 , tf , N ); dt = t( 2 ) - t( 1 );
% Prealocarea si valori de inceput
vx = zeros( 1 , N ); vy = vx ; % prealocare pt componentele vitezei 
x = zeros( 1 , N ); y = x ; % prealocarea pt coordonatele distantei 
vx( 1 ) = v0 * cosd( alpha0 );
vy( 1 ) = v0 * sind( alpha0 );

% simularea
for i = 1:N-1
    aux = 1 - dt *( r1 + r2 * sqrt( vx( i )^2 + vy( i )^2) ) / m ;
    vx( i + 1 ) = vx( i ) * aux ; % recurenta de ordin 2 pt vx
    vy( i + 1 ) = vy( i ) * aux  - g*dt; % recurenta de ordin 2 pt vy
    
    x( i + 1 ) = x( i ) + vx( i ) * dt ; % miscare uniforma pe intervalul dt 
    y( i + 1 ) = y( i ) + vy( i ) * dt ; % miscare uniforma pe intervalul dt
    
    if y( i + 1 ) < 0 
        break;
    endif 
endfor 

%eliminarea zerourilor de la final
t = t( 1:i );
vx = vx( 1:i );
vy = vy( 1:i );
x = x( 1:i );
y = y( 1:i );

% reprezentarea grafica 
if Gv == 1 
  figure( 1 );
  plot( t , vx , '-r' , t , vy , '-b');
  xlabel( 't[s]' ); ylabel( 'v[m/s]' );
  grid;
  title('Componentele vitezei ca functii de timp');
  legend( 'Vx' , 'Vy' );
endif  
  
if Gc == 1 
  figure(2);
  plot( t , x/1e3 , '-r' , t , y/1e3 , '-b');
  xlabel( 't[s]' ); ylabel( 'coordonata[km]' );
  grid;
  title('Coordonatele ca functii de timp');
  legend( 'x' , 'y' , 'Location' , 'NorthWest');
endif 

if T == 1 
  figure(3);
  plot( x/1e3 , y/1e3 , '-k' , 'LineWidth' , 2 );
  xlabel( 'x[km]' ); ylabel( 'y[km]' );
  grid;
  title('Curba balistica');
   
endif 

% afisarea unei marimi de interes 
tf = t( i ); % timpul de zbor 
b = x( i ); % bataia
h = max( y ); % altitudinea maxima
tu = t( y == h ); % timpul de urcare
tc = tf - tu ;  
Q = m * 0.5 * ( v0^2 -  vx( i )^2 - vy( i )^2 ); % Caldura disipata
afis=[ 'Timpul de zbor: ', num2str( tf ) , ' s']; disp( afis );
afis=[ 'Bataia proiectilului: ', num2str( b/1e3 ) , ' km']; disp( afis );
afis=[ 'Altitudinea maxima: ', num2str( h/1e3 ) , ' km']; disp( afis );
afis=[ 'Timpul de urcare: ', num2str( tu ) , ' s']; disp( afis );
afis=[ 'Timpul de coborare: ', num2str( tc ) , ' s']; disp( afis );
afis=[ 'Caldura produsa: ', num2str( Q/1e3 ) , ' kJ']; disp( afis );
 
 if F == 1 
    figure(4);
    for j = 1:i 
      plot( x/1e3 , y/1e3 , '-c' ); hold on;
      xlabel( 'x[km]' ); ylabel( 'y[km]' );
      grid;
      title('Animatie curba balistica');
      axis equal; axis tight;
      plot( x( j )/1e3 , y( j )/1e3 , 'ok'  );
    endfor
 endif
