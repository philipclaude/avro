clear;clc;close all;

d = 2;
fid = fopen(['../../tmp/bluenoise-dim',num2str(d),'-n1000.json']); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
data = jsondecode(str);
X = data.x;

n = size(X,1)/d;
x = X(1:d:end);
y = X(2:d:end);
z = X(3:d:end);
t = X(4:d:end);

X = ceil(n*x);
Y = ceil(n*y);
Z = int32(n*z);
T = int32(n*t);

D = zeros(n,n);

for i=1:n
    D(X(i),Y(i)) = 1;
end

f  = fft2(D);
fs = fftshift(f);
ps = log(abs(fs));

figure;
imshow(D);

figure;
imshow(ps);

% compute the rdf and power spectrum
fidtxt = fopen('points.txt','w');
fprintf(fidtxt,'%d\n',n);
for i=1:n
    fprintf(fidtxt,'%g %g\n',x(i),y(i));
end
fclose(fidtxt);
system('~/Codes/psa/psa --rp points.txt');
%system('pdflatex points_rp.tex');
