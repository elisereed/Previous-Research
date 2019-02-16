function [ftilda, x] = RP(x, N, f, sigma, mu, noise)

%line search
opts = optimset('Display', 'off', 'LargeScale', 'off', 'TolX', mu);
funeval = noise(sigma);
n = size(x);

k = 1;

while k < N
    
    u = randn(n);
    u = u/norm(u);
    %uniformly random vector from S_n-1
    
    addnoise = noise(sigma);
    
    noisyf =@(x) f(x) + addnoise;
    
    [h, ftilda, ~, infos] = ...
        fminunc(@(h) feval(noisyf, x + h*u), 0, opts);
    
    funeval = funeval + infos.funcCount;
    x = x + h*u;
    
    k = k+1;
end

    
    
    