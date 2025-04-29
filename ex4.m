%% ex4 
% Solving ODE (k(x)*u(x)')' + c(x)u(x) = f(x), u(0) = u(1) = 0 using CG 
% with respect to H_0^1 inner product, then compare with pcg
% where f(x), c(x) are piecwise continuous 
clear all

dom = [0 0.5 1]; 

% Define c(x) as a piecewise function
c = chebfun({@(x) 1, @(x) 10}, dom);

% Plot c(x)
figure(1)
plot(c, 'LineWidth', 2); 
grid on
title('c(x)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('x', 'FontSize', 11);
box on

% Export the figure to a PDF file
exportgraphics(gcf, 'discont_c.pdf', 'ContentType', 'vector', ...
    'BackgroundColor', 'none')

% Define f(x) as a piecewise function
f = chebfun({@(x) -6*x + x.^3, @(x) -30*x.^3 + 40*x.^2 + 8*x - 8}, dom);

% Plot f(x)
figure(2)
plot(f, 'LineWidth', 2);
grid on
title('f(x)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('x', 'FontSize', 11);
box on

% Export the figure to a PDF file
exportgraphics(gcf, 'discont_f.pdf', 'ContentType', 'vector', ...
    'BackgroundColor', 'none')

x = chebfun('x', [0,1]);
k = chebfun(0*x + 1);
u0 = chebfun(0*x);
iter = 7;
exact = chebfun({@(x) x.^3, @(x) -3*x.^3 + 4*x.^2 - x}, dom);

% Plot exact solution
figure(3)
plot(exact, 'LineWidth', 2); 
grid on
title('u(x)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('x', 'FontSize', 11);
box on

% Export the figure to a PDF file
exportgraphics(gcf, 'discont_u.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')

% Compute aproximations using pcg and cgh1_ab

[u,sols1] = cgh1_ab(k,c,f,u0,iter,0,1);

L = chebop(@(u) -diff(k .* diff(u)) + c .* u,dom);
[u2,flag,relres,iter,resvec2,sols2] = pcg2(L,f);


errors1 = zeros(length(sols1),1);
for i = 1:length(sols1)
    err = sols1{i}-exact;
    errors1(i) = sqrt(sum(k*diff(err)*diff(err) + c*(err)*(err)));
end

errors2 = zeros(length(sols2),1);
for i = 1:length(sols2)
    err = sols2{i}-exact;
    errors2(i) = sqrt(sum(k*diff(err)*diff(err) + c*(err)*(err)));
end

% Plot errors 
figure(4)
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
exportgraphics(gcf, 'discont_err.pdf', 'ContentType', 'vector', ...
    'BackgroundColor', 'none')
