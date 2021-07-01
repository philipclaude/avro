function metric
clear;clc;close all;

% case 0a: my favourite metric in [-5,5]^2
[mx,xlims,mr,rlims] = sshock2d;
V0a = volume( sqrt(det(mx)),xlims);
fprintf('sshock case expects %d triangles in [%d,%d]\n',V0a/volume0(2),V0a/(volume0(2,sqrt(2.))),V0a/volume0(2,1./sqrt(2.)));

% case 0b: my favourite metric in C([5,5],5)
s = symvar(mr);
%assert( s(1)=='r' && s(2)=='theta' );
V0b = volume( sqrt(det(mr))*s(1) ,rlims);
fprintf('sshock case expects %d triangles in [%d,%d]\n',V0b/volume0(2),V0b/(volume0(2,sqrt(2.))),V0b/volume0(2,1./sqrt(2.)));

% case 0c: boundarylayer2d
[mx,xlims] = boundarylayer2d;
V0c = volume( sqrt(det(mx)),xlims);
fprintf('bl2d case expects %d triangles in [%d,%d]\n',V0c/volume0(2),V0c/(volume0(2,sqrt(2.))),V0c/volume0(2,1./sqrt(2.)));

% case 1a: linear metric in box [0,1]^3
[m,lims] = linearcube3d;
V1 = volume(sqrt(det(m)),lims)*1*1; % pedantic, we only integrated in z
fprintf('CL = %g\n',V1);
fprintf('linear cube case expects %d tetrahedra in [%d,%d]\n',V1/volume0(3),V1/(volume0(3,sqrt(2.))),V1/volume0(3,1./sqrt(2.)));

% case 1b: linear metric in tesseract [0,1]^4
[m,lims] = linearcube4d;
V1b = volume(sqrt(det(m)),lims)*1*1; % pedantic, we only integrated in z
fprintf('TL = %g\n',V1b);
fprintf('linear tesseract case expects %d pentatopes in [%d,%d]\n',V1b/volume0(4),V1b/(volume0(4,sqrt(2.))),V1b/volume0(4,1./sqrt(2.)));

% case 2: linear metric in box [0,1]^3 -C([0,0,0],0.5) (cylinder)
% since there is no variation in x or y, the number of elements
% is the ratio of the volume of (box -cyl)/box
Vbox = 1.; % unit cube
Vcyl = 0.25*pi*0.5^2*1; % quarter cylinder, r = 0.5, l = 1.
V2 = V1*( Vbox -Vcyl )/Vbox; % volume ratio
fprintf('linear cube-cyl case expects %d tetrahedra in [%d,%d]\n',V2/volume0(3),V2/(volume0(3,sqrt(2.))),V2/volume0(3,1./sqrt(2.)));

% case 3: polar 1
[mx,xlims,mr,rlims] = polar1;
s = symvar(mr);
assert(s(1)=='r' && s(2)=='t');
r = s(1);
Vr = volume( sqrt(det(mr))*r , rlims )*1; % integration in r,t
Vx = volume( sqrt(det(mx)) , xlims )*1; % integration in x,y over box
V3 = Vx -Vr;
fprintf('polar1 cube case expects %d tetrahedra in [%d,%d]\n',Vx/volume0(3),Vx/(volume0(3,sqrt(2.))),Vx/volume0(3,1./sqrt(2.)));
fprintf('polar1 cube-cyl case expects %d tetrahedra in [%d,%d]\n',V3/volume0(3),V3/(volume0(3,sqrt(2.))),V3/volume0(3,1./sqrt(2.)));

% case 4: polar 2, integrate to within +/- 10 (doesn't matter that much)
[mx,xlims,mr,rlims] = polar2;
s = symvar(mr);
assert(s(1)=='r' && s(2)=='theta');
r = s(1);
Vr = volume( sqrt(det(mr))*r , rlims , 10 )*1; % integration in r,t
Vx = volume( sqrt(det(mx)) , xlims , 10  )*1; % integration in x,y over box
V4 = Vx -Vr;
fprintf('polar2 cube case expects %d tetrahedra in [%d,%d]\n',Vx/volume0(3),Vx/(volume0(3,sqrt(2.))),Vx/volume0(3,1./sqrt(2.)));
fprintf('polar2 cube-cyl case expects %d tetrahedra in [%d,%d]\n',V4/volume0(3),V4/(volume0(3,sqrt(2.))),V4/volume0(3,1./sqrt(2.)));

end

function v0 = volume0(d,a)
% volume of unit simplex
if nargin==1
    a = 1.;
end
v0 = a^d*sqrt(d+1.)/(factorial(d)*sqrt(2.^d));
end

function V = volume(dV,lims,tol)

% +/- 10ish elements who really cares
if nargin==2
    tol = 1;
end

% integrate in each direction
nd = length(lims);
if nd==1
    bounds1 = lims{1};
    fh = matlabFunction(dV);
    V = integral(fh,bounds1(1),bounds1(2),'abstol',tol);
elseif nd==2
    bounds1 = lims{1};
    bounds2 = lims{2};
    fh = matlabFunction(dV);
    V = integral2(fh,bounds1(1),bounds1(2),bounds2(1),bounds2(2),'abstol',tol);
    %V = dblquad(fh,bounds1(1),bounds1(2),bounds2(1),bounds2(2),'abstol',tol);
elseif nd==3
    bounds1 = lims{1};
    bounds2 = lims{2};
    bounds3 = lims{3};
    fh = matlabFunction(dV);
    V = integral3(fh,bounds1(1),bounds1(2),bounds2(1),bounds2(2),bounds3(1),bounds3(2),'abstol',tol);
end

end

function [mx,xlims] = boundarylayer2d
y = sym('y','real');
hy = 1e-4 +0.1*y;
mx(2,2) = 1./(hy*hy);
mx(1,1) = 1;
mx(1,2) = 0;
mx(2,1) = 0;

xlims = {[0,1]};
end

function [mx,xlims,mr,rlims] = sshock2d

x = sym('x','real');
y = sym('y','real');

H  = 5.;
f  = 10.;
A  = 2.;
P  = 5.;
y0 = 5.;
w  = pi/P;

phi  = -f*( (y+5) -A*cos(w*(x+5)) -y0 );
z    = H*tanh( -phi );
dzdx = diff(z,x);
dzdy = diff(z,y);

mx(1,1) = 1. +dzdx*dzdx;
mx(1,2) = dzdx*dzdy;
mx(2,1) = mx(1,2);
mx(2,2) = 1. +dzdy*dzdy;

xlims = {[-15,15],[-5,5]};

r = sym('r','real');
theta = sym('theta','real');

mr = mx;
rlims = {[0,20],[0,2*pi]};

end

function [m,lims,vars] = linearcube3d
hx = 0.1;
hy = 0.1;
h0 = 1e-3;
z  = sym('z','real');
hz = h0 +2*(0.1 -h0)*abs(z-0.5);
m = [1./hx^2,0,0;0,1/hy^2,0;0,0,1/hz^2];
lims = {[0,1]};
vars = z;
end

function [m,lims,vars] = linearcube4d
hmin = 0.0025;
hx = 100*hmin;
hy = 100*hmin;
hz = 100*hmin;
h0 = hmin;
t  = sym('t','real');
ht = h0 +2*(hx -h0)*abs(t-0.5);
m = [1./hx^2,0,0,0;0,1/hy^2,0,0;0,0,1/hz^2,0;0,0,0,1./ht^2];
lims = {[0,1]};
vars = t;
end

function [mx,xlims,mr,rlims] = polar1
r = sym('r','real');
t = sym('t','real');
ht = 0.1;
hz = 0.1;
h0 = 1e-3;
hr = h0 +2*(0.1-h0)*abs( r-0.5 );

Q = [cos(t),sin(t),0;-sin(t),cos(t),0;0,0,1];
D = [1./hr^2,0,0;0,1/ht^2,0;0,0,1/hz^2];
mr = Q'*D*Q;
rlims = {[0.,.5],[0,.5*pi]};

x = sym('x','real');
y = sym('y','real');
r = sqrt( x*x +y*y );
cost = x/r;
sint = y/r;
Q = [cost,sint,0;-sint,cost,0;0,0,1];
hr = h0 +2*(0.1 -h0)*abs( r -0.5 );
D = [1./hr^2,0,0;0,1/ht^2,0;0,0,1/hz^2];
mx = Q'*D*Q;
xlims = {[0,1],[0,1]};

f = 2;
mr = f*mr;
mx = f*mx;
end

function [mx,xlims,mr,rlims] = polar2
H = 100.; % it doesn't matter what this is because it gets subtracted off
r = sym('r','real');
theta = sym('theta','real');
hz = 0.1;
h0 = 1e-3;
d = 10*(0.6 -r);
ht = .5*(1-sign(d))*0.1 +.5*(1+sign(d))*(d/40. +0.1*(1.-d));
R = r -0.5;
ht = .5*(1+sign(R))*ht +H*(1-sign(R)); % will be 2*H when r < 0.5, but ht otherwise

hr = h0 +2*(0.1-h0)*abs( r-0.5 );

Q = [cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1];
D = [1./hr^2,0,0;0,1/ht^2,0;0,0,1/hz^2];
mr = Q'*D*Q;
rlims = {[0.,.5],[0,.5*pi]};

x = sym('x','real');
y = sym('y','real');
r = sqrt( x*x +y*y );
cost = x/r;
sint = y/r;
Q = [cost,sint,0;-sint,cost,0;0,0,1];

%d = 10*(0.6 -r);
%ht = .5*(1-sign(d))*0.1 +.5*(1.+sign(d))*(d/40. +0.1*(1.-d));
%R = r -0.5;
%ht = .5*(1+sign(R))*ht +H*(1-sign(R)); % will be 2*H when r < 0.5, but ht otherwise

R  = r - 0.5;
d0 = 10*abs(r-0.5);
ht = 0.1*d0 + 0.025*(1.0 - d0 );

hr = h0 +2*(0.1 -h0)*abs( r-0.5 );
D = [1./hr^2,0,0;0,1/ht^2,0;0,0,1/hz^2];
mx = Q'*D*Q;
xlims = {[0,1],[0,1]};

end
