function [ftilda, x] = STARS(x, N, f, sigma, L1,noise)


n = length(x);
ftildaold = feval(f,x) + feval(noise,sigma);

%created separate script for noise function

mu = ((8*sigma^2*n)/(L1^2*(n+6)^3))^.25;
%constant additive smoothing stepsize

h = 1/(4*(n+4)*L1);
%constant additive stepsize
k=1;

while k < N
    
      
    u = normrnd(0,1,n,1);
    %nx1 vector of iid normally distributed random numbers
   
    ftildanew = feval(f,x + mu*u) + feval(noise,sigma);
    
    s = ((ftildanew-ftildaold)/mu).*u;
    
    x = x - h*s;
    
    ftildaold = feval(f,x) + feval(noise,sigma);
    
    k = k+1;
    
end

ftilda = ftildaold;
return
    