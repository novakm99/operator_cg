%% ex3
% Compare absolute errors of functions cgh1_ab and pcg2 
% for a specific example of ODE 
% (k(x)*u(x)')' + c(x)u(x) = f(x) on (a,b), u(a) = 1, u(b) = 3
%
% Solving the auxiliary problem (k(x)*v(x)')' + c(x)v(x) = modified f(x),
% v(0) = v(1) = 0 using Conjugate Gradient (CG) with respect to 
% the H_0^1 inner product, then recovering u(x) = v(x) + z(x) to satisfy
% original boundary conditions and compare with pcg2.
clear all

% Define the variable x over the interval [-1, 1] using Chebfun
x = chebfun('x', [0,1]);

% Initialization
iter =16;
f = chebfun(x.^4 + x.^3 + x.^2 - 4*x -3);
k = chebfun(1 + x);
c = chebfun(x.^2);
u0 = chebfun(0, [0,1]);
z = chebfun(1 + 2*x);
[u,sols1] = cgh1_ab(k,c,f + 2 - c*z,u0,iter,0,1);

exact = chebfun(x.^2 + x +1);

L = chebop(@(u) -diff(k .* diff(u)) + c .* u, [0,1]);
L.lbc = 1;
L.rbc = 3;

[u2,~,~,~,~,sols2] = pcg2(L,f);

% Calculate absolute errors for each solution from cgh1_ab and pcg2
errors1 = zeros(length(sols1),1);
for i = 1:length(sols1)
    err = sols1{i} + z - exact;
    errors1(i) = sqrt(sum(k*diff(err)*diff(err)+ c*(err)*(err)));
end

errors2 = zeros(length(sols2),1);
for i = 1:length(sols2)
    err = sols2{i}-exact;
    errors2(i) = sqrt(sum(k*diff(err)*diff(err) + c*(err)*(err)));
end

% Plot errors
figure()
semilogy(errors1, '-', 'LineWidth', 1.8)
hold on
semilogy(errors2,'-', 'LineWidth', 1.8)
legend('cgh1','pcg')
xlabel('Iterace', 'FontSize', 14);
hold off
grid on;
set(gca, 'FontSize', 12);
box on;

% Export the figure to a PDF file
exportgraphics(gcf, 'H01cg_err_ab.pdf', 'ContentType', 'vector', ...
    'BackgroundColor', 'none')