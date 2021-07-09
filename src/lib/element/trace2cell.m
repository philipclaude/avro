clear;clc;close all;

nCell  = 4; % nCell-simplex
nTrace = nCell -1;

v = 0:nTrace;
A = perms(v);
A = flipud(A); % reverse the number of the permutations

np = size(A,1);
d  = zeros(np,1);

s0 = sym('s0','real');
t0 = sym('t0','real');
u0 = sym('u0','real');

if nTrace==1 % line trace
    S0 = [1-s0,s0]';
elseif nTrace==2 % triangle trace
    S0 = [1-s0-t0,s0,t0]';
elseif nTrace==3 % tetrahedron trace
    S0 = [1-s0-t0-u0,s0,t0,u0]';
end

Y = eye(nCell);
y = zeros(nTrace,1);

% reference coordinates
M = [zeros(nTrace,1),eye(nTrace)];

for i=1:np
    
    % extract the permutation
    x = A(i,:);
        
    % compute the levi-cevita sign (orientation)
    B = Y(:,x+1);
    d(i) = det( B );
    
    if d(i)<0
        orientation = 'negative';
    else
        orientation = 'positive';
    end
    fprintf('permutation %d (%s):',i,orientation);
    disp(x);
    
    % compute the mapping from original barycentrics to permuted simplex
    y(x+1) = v;
    
    % display the coordinates
    S = M(:,y+1)*S0;
    disp(S);
 
end

A = [A,d];
disp(A);
