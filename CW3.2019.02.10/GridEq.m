function z = GridEq(n,a,b)

% z=GridEq(n,a,b) computes n+1 equally spaced points in [a,b]. (z(1)=a;
% z(n+1)=b)

% Initialize the vecter z
z = zeros(1,n+1);

for j = 0:n
    z(j+1) = a+(b-a)*j/n;
end
    
    
    
