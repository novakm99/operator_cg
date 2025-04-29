function [u,sols] = cgl2(k,c,f,u0,maxit,dim)
% The Conjugate gradient method using theory of Riesz map with respect to 
% L^2 inner product.
%
% Solving ODE:
% (k(x)*u(x)')' + c(x)u(x) = f(x), on (-1,1); u(-1)=u(1)=0
%
% Inputs:
% k > 0 on (-1,1)
% c >= 0 on (-1,1)
% f  
% u0  - initial guess for the solution, must be from H_0^1
% maxit - maximum number of iterations
% dim - dimension of the projection space consisting of shifted Chebyshev 
% polynomials that satisfy the boundary conditions
%
% Outputs:
% u - final solution after the iterations
% sols - a cell array containing the solutions at each iteration
% NOTE: This code requires the chebfun package to handle Chebyshev 
% functions and operations over intervals. If chebfun is not installed,
% you can download it from https://www.chebfun.org/ 

sols = cell(maxit+1,1);  
sols{1} = u0; 

tau_b = myprojection(f,dim); %do projection
tau_Au = myprojection(-diff(k*diff(u0))+c*u0,dim); %do projection
tau_r = tau_b - tau_Au;
p = tau_r;
u = u0;
for n =1:maxit
    tau_Ap = myprojection(-diff(k*diff(p))+c*p,dim); %do projection
    alpha = sum(tau_r*tau_r)/sum(tau_Ap*p);
    u = u + alpha * p;
    sols{n+1} = u; 
    q = sum(tau_r*tau_r);
    tau_r = tau_r - alpha*tau_Ap;
    beta = sum(tau_r*tau_r)/q;
    p = tau_r + beta*p;
end