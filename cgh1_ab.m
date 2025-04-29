function [u,sols] = cgh1_ab(k,c,f,u0,maxit,a,b)
% The Conjugate gradient method using theory of Riesz map with respect to 
% H_0^1 inner product.
%
% Generalization of function cgh1 for solving the boundary value problem:
% (k(x)*u(x)')' + c(x)u(x) = f(x), on (a,b); u(a)=u(b)=0
%
% Inputs:
% k > 0 
% c >= 0
% f  
% u0  - initial guess for the solution, must be from H_0^1
% maxit - maximum number of iterations
% k,c,f,u0 must be Chebfuns of variable 'x' defined on  [a,b]
%
% Outputs:
% u - final solution after the iterations
% sols - a cell array containing the solutions at each iteration
% NOTE: This code requires the chebfun package to handle Chebyshev 
% functions and operations over intervals. If chebfun is not installed,
% you can download it from https://www.chebfun.org/ 

x = chebfun('x', [a,b]);

% Cell array to store the solutions at each iteration
sols = cell(maxit+1,1);  
sols{1} = u0;

tau_b = -cumsum(cumsum(f)) + sum(cumsum(f))*chebfun(x-a)/(b-a);
g = cumsum(k*diff(u0));
h = cumsum(cumsum(c * u0));
tau_Au = g - 1/(b-a)*g(b)*chebfun(x-a) - h + 1/(b-a)*h(b)*chebfun(x-a);
tau_r = tau_b - tau_Au;
p = tau_r;
u = u0;

% Iterative process: Conjugate Gradient method
for n =1:maxit
    g = cumsum(k * diff(p));
    h = cumsum(cumsum(c*p));
    tau_Ap = g - 1/(b-a)*g(b)*chebfun(x-a) - h + 1/(b-a)*h(b)*chebfun(x-a);
    sum(diff(tau_Ap)*diff(p));
    alpha = sum(diff(tau_r)*diff(tau_r))/sum(diff(tau_Ap)*diff(p));
    u = u + alpha * p;
    sols{n+1} = u;
    q = sum(diff(tau_r)*diff(tau_r));
    tau_r = tau_r - alpha*tau_Ap;
    beta = sum(diff(tau_r)*diff(tau_r))/q;
    p = tau_r + beta*p;
end