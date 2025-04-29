%% ex1
% Compare residuals, absolute errors and numbers of Chebyshev points
% of functions cgh1 and pcg2 for solving a specific example of ODE 
% (k(x)*u(x)')' + c(x)u(x) = f(x), u(-1) = u(1) = 0

clear all

% Define variable x over the interval [-1, 1] using chebfun
x = chebfun('x', [-1,1]);

% Initialization
maxit = 23;
f = chebfun(-pi*cos(pi*x) + pi.^2*(x + 3)*sin(pi*x));
k = chebfun(x + 2);
c = chebfun(pi^2);
u0 = chebfun(0);


% Define exact solution for comparison
exact = chebfun(sin(pi*x));

[u1,sols1] = cgh1(k,c,f,u0,maxit);


L = chebop(@(u) -diff(k .* diff(u)) + c .* u);
[u2,~,~,~,~,sols2] = pcg2(L,f);


% Calculate absolute errors for each solution from cgh1 and pcg2
errors1 = zeros(length(sols1),1);
for i = 1:length(sols1)
    err = sols1{i}-exact;
    errors1(i) = sqrt(sum(k*diff(err)*diff(err)+ c*(err)*(err)));
end

errors2 = zeros(length(sols2),1);
for i = 1:length(sols2)
    err = sols2{i}-exact;
    errors2(i) = sqrt(sum(k*diff(err)*diff(err)+ c*(err)*(err)));
end

% Calculate residual norms for cgh1 and pcg2 solutions
resvec1 = zeros(length(sols1),1);
for i = 1:length(sols1)
    resvec1(i) = norm(L(sols1{i})-f ,2);
end

resvec2 = zeros(length(sols2),1);
for i = 1:length(sols2)
    resvec2(i) = norm(L(sols2{i})-f ,2);
end

% Computing numbers of Chebyshev points
dims1 = zeros(length(sols1),1);
for i = 1:length(sols1)
    dims1(i) = length(sols1{i});
end

dims2 = zeros(length(sols2),1);
for i = 1:length(sols2)
    dims2(i) = length(sols2{i});
end

% Plot the residuals
figure(1)
semilogy(resvec1,'-', 'LineWidth', 1.8)
hold on
semilogy(resvec2,'-', 'LineWidth', 1.8)
legend('cgh1','pcg')
xlabel('Iterace', 'FontSize', 14);
hold off
grid on;
set(gca, 'FontSize', 12);
box on;
% Export the figure to a PDF file
exportgraphics(gcf, 'H01cg_res.pdf', 'ContentType', 'vector', ...
    'BackgroundColor', 'none')
%%
% Plot the errors 
figure(2)
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
exportgraphics(gcf, 'H01cg_err.pdf', 'ContentType', 'vector', ...
    'BackgroundColor', 'none')
%%

% Plot numebrs of Chebyshev points 
figure(3)
plot(dims1,'.-', 'LineWidth', 2, 'MarkerSize', 18)
hold on
plot(dims2,'.-', 'LineWidth', 2, 'MarkerSize', 18)
legend('cgh1','pcg', 'FontSize', 12, 'Location', 'southeast')
grid on;
set(gca, 'FontSize', 12);
xlabel('Iterace', 'FontSize', 14);
box on;
hold off

% Export the figure to a PDF file
exportgraphics(gcf, 'H01cg_dim.pdf', 'ContentType', 'vector', ...
    'BackgroundColor', 'none')
