function [ftilda, x] = RG(x, N, f, sigma, L1, noise)

n = length(x);
h = 1/(4*(n+4)*L1);
mu = (5/(3*(n+4)))*sqrt((2^-16)/(2*L1));

k =1;

while k < N
    
    u = normrnd(0,1,n,1);
    %nx1 vector of normally distributed random numbers
    
    ftildaold = feval(f,x) + feval(noise,sigma);
    ftildanew = feval(f,x+mu*u) + feval(noise,sigma);
    
    s = ((ftildanew-ftildaold)/mu)*u;
    
    x = x - h*s;
    
    k = k+1;
end
ftilda = ftildanew;
return