clear;clc;close all;

d = 2;
fid = fopen(['bluenoise-dim',num2str(d),'-n1000.json']); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
data = jsondecode(str);

r = data.r;
X = data.x;
n = size(X,1)/d;
x = X(1:d:end);
y = X(2:d:end);
z = X(3:d:end);
t = X(4:d:end);

X = int32(n*x);
Y = int32(n*y);
Z = int32(n*z);
T = int32(n*t);

if d == 2
    D = zeros(n,n);

    for i=1:n
        D(X(i),Y(i)) = 1;
    end
    
    for i=1:n
        for j=1:n
            D(i,j) = (x(i) - x(j))^2 + (y(i)-y(j))^2;
        end
    end
    
elseif d == 3
    D = zeros(n,n,n);
    for i=1:n
        D(X(i),Y(i),Z(i)) = 1;
    end
elseif d == 4
    D = zeros(n,n);
    for i=1:n
        D(X(i),Y(i)) = 1;
    end
end

f  = fftn(D);
fs = fftshift(f);
ps = log(abs(fs));

figure;
imshow(D);

figure;
imshow(ps);

return;

% compute the power spectral density
nannuli = 10;

var = zeros(nannuli,1);
step = 1.0/nannuli;
domain2 = int32(n/2.0);

power = zeros(nannuli,1);
anisotropy = zeros(nannuli,1);

for i=1:nannuli
    
    rmin = int32( (sqrt(2.0)*domain2 * i     * step)^2 );
    rmax = int32( (sqrt(2.0)*domain2 * (i+1) * step)^2 );
    
    Nr = 0;
    for y=1:n
        for x=1:n
            r = int32( (x-domain2)^2 + (y-domain2)^2 );
            if (r >= rmin && r <= rmax)
                power(i) = power(i) + ps(x,y);
                Nr = Nr +1;
            end
        end
    end
end

plot(power);