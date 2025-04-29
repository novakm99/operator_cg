function [u,sols] = cgh1(k,c,f,u0,maxit)
% The Conjugate gradient method using theory of Riesz map with respect to 
% H_0^1 inner product.
%
% Solving ODE:
% (k(x)*u(x)')' + c(x)u(x) = f(x), on (-1,1); u(-1)=u(1)=0
%
% Inputs:
% k > 0 
% c >= 0
% f  
% u0  - initial guess for the solution, must be from H_0^1
% maxit - maximum number of iterations
% k,c,f,u0 must be Chebfuns of variable 'x' defined on (-1,1)
%
% Outputs:
% u - final solution after the iterations
% sols - a cell array containing the solutions at each iteration
% NOTE: This code requires the chebfun package to handle Chebyshev 
% functions and operations over intervals. If chebfun is not installed,
% you can download it from https://www.chebfun.org/ 

% Cell array to store the solutions at each iteration
sols = cell(maxit+1,1);  
sols{1} = u0;  

% Define 'x' as a Chebyshev function over the interval [-1, 1]
x = chebfun('x', [-1,1]); 

tau_b = -(cumsum(cumsum(f)) - 0.5*sum(cumsum(f))*chebfun(x+1));
g = cumsum(k*diff(u0));  
h = cumsum(cumsum(c * u0)); 
tau_Au = g - 0.5*g(1)*chebfun(x+1) - h + 0.5*h(1)*chebfun(x+1);
tau_r = tau_b - tau_Au;
p = tau_r;
u = u0;

% Iterative process: Conjugate Gradient method
for n = 1:maxit
    g = cumsum(k * diff(p));  
    h = cumsum(cumsum(c * p));  
    tau_Ap = g - 0.5*g(1)*chebfun(x+1) - h + 0.5*h(1)*chebfun(x+1);
    alpha = sum(diff(tau_r)*diff(tau_r)) / sum(diff(tau_Ap)*diff(p));
    u = u + alpha * p;
    sols{n+1} = u;
    q = sum(diff(tau_r)*diff(tau_r)); 
    tau_r = tau_r - alpha * tau_Ap; 
    beta = sum(diff(tau_r)*diff(tau_r)) / q;
    p = tau_r + beta * p;
end