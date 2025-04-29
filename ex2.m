%% ex2
% Computing approximations of the solution to the specific ODE:
% (k(x) * u'(x))' + c(x) * u(x) = f(x), with boundary conditions 
% u(-1) = u(1) = 0.
% This is done using Conjugate Gradient (CG) method with respect to the L^2 
% inner product, where projections are applied to solve the problem 
% in finite-dimensional subspaces.
%
% The goal is to compare the absolute errors for different dimensions of 
% the subspace, where the solution is projected (e.g., dim = 5, 10, 20, 30)

clear all

% Define the variable x over the interval [-1, 1] using Chebfun
x = chebfun('x', [-1,1]);

% Initialization
maxit = 60;

f = chebfun(-pi*cos(pi*x) + pi.^2*(x + 3)*sin(pi*x));
k = chebfun(x + 2);
c = chebfun(pi^2);
u0 = chebfun(0);

[u5,sols5] = cgl2(k,c,f,u0,maxit,5);
[u10,sols10] = cgl2(k,c,f,u0,maxit,10);
[u20,sols20] = cgl2(k,c,f,u0,maxit,20);
[u30,sols30] = cgl2(k,c,f,u0,maxit,30);

% Define exact solution
exact =  chebfun(sin(pi*x));

% Computing errors for each case (5, 10, 20, 30 dimensions)
errors5 = zeros(length(sols5),1);
for i = 1:length(sols5)
    err = sols5{i}-exact;
    errors5(i) = sqrt(sum(k*diff(err)*diff(err) + c*(err)*(err)));
end

errors10 = zeros(length(sols10),1);
for i = 1:length(sols10)
    err = sols10{i}-exact;
    errors10(i) = sqrt(sum(k*diff(err)*diff(err) + c*(err)*(err)));
end

errors20 = zeros(length(sols20),1);
for i = 1:length(sols20)
    err = sols20{i}-exact;
    errors20(i) = sqrt(sum(k*diff(err)*diff(err) + c*(err)*(err)));
end

errors30 = zeros(length(sols30),1);
for i = 1:length(sols30)
    err = sols30{i}-exact;
    errors30(i) = sqrt(sum(k*diff(err)*diff(err) + c*(err)*(err)));
end

% Plot errors
figure
semilogy(errors5,'LineWidth', 1.8)
hold on
semilogy(errors10,'LineWidth', 1.8)
semilogy(errors20,'LineWidth', 1.8)
semilogy(errors30,'LineWidth', 1.8)
legend('dim=5','dim=10','dim=20','dim=30', 'Location', 'southwest')
xlabel('Iterace', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 12);
box on;

% Export the figure to a PDF file
exportgraphics(gcf, 'L2cg2.pdf', 'ContentType', 'vector', ...
    'BackgroundColor', 'none')