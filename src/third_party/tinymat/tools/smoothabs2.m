syms EPS

b = [ 0;EPS/2;0;1;0;0];
 
A = [  1, (-EPS/2)^1,   (-EPS/2)^2,   (-EPS/2)^3,    (-EPS/2)^4,    (-EPS/2)^5 ;
       1, ( EPS/2)^1,   ( EPS/2)^2,   ( EPS/2)^3,    ( EPS/2)^4,    ( EPS/2)^5 ;
       0,          1, 2*(-EPS/2)^1, 3*(-EPS/2)^2,  4*(-EPS/2)^3,  5*(-EPS/2)^4 ;
       0,          1, 2*( EPS/2)^1, 3*( EPS/2)^2,  4*( EPS/2)^3,  5*( EPS/2)^4 ;
       0,          0,            2, 6*(-EPS/2)^1, 12*(-EPS/2)^2, 20*(-EPS/2)^3 ;
       0,          0,            2, 6*( EPS/2)^1 12*( EPS/2)^2, 20*( EPS/2)^3 ];
  
     
alpha = A\b;

syms X

alphaT = A(1:4,1:4)\b(1:4);

% function with coefficients
Alpha = matlabFunction(flip(alpha));
AlphaT = matlabFunction(flip(alphaT)); 

clear EPS
y = @(x,EPS) polyval(Alpha(EPS),x);
yT = @(x,EPS) polyval(AlphaT(EPS),x);

EPS = 0.5;
X = linspace(-EPS,EPS);
clf
plot(X,y(X,EPS),'r') % quintic
hold on
plot(X,yT(X,EPS),'k:') % cubic
plot(X,max(0,X),'b--')