function z = GridCheb(n,a,b)

% z=GridCheb(n,a,b) computes n+1 Chebychev points in [a,b]

% Initialize the vecter z
z = zeros(1,n+1);

for j = 0:n
    z(j+1) = 1/2*(a+b)+1/2*(b-a)*cos((2*j+1)/(2*n+2)*pi);
end