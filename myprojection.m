function proj = myprojection(f, dim)
% Projects a function f onto the span of shifted Chebyshev 
% polynomials such that they are zero in -1, 1
%   Inputs:
%       f   - chebfun object (function to project)
%       dim - number of basis functions to use
%   Output:
%       proj - the Chebyshev projection of f

% Ensure the domain matches f's domain
dom = domain(f);

% Generate the first 'dim' even-degree Chebyshev polynomials, shifted by -1
T = cell(dim, 1);

num_even = floor(dim / 2);
num_odd = dim - num_even;
for k = 1:num_even
    T{k} = chebpoly(2*k, dom) - 1;  % Shifted even Chebyshev polynomial
    T{k} = T{k}/norm(T{k});
end
for k =1:num_odd
    T{num_even + k} = chebpoly(2*k+1,dom) -chebfun(@(x) x); % Shifted odd Chebyshev polynomial
end

% Calculate Gram matrix G where G(i,j) = (P_i, P_j)_{L^2}
    G = zeros(dim, dim);
    for i = 1:dim
        for j = 1:dim
            G(i,j) = sum(T{i} .* T{j});  % L^2 inner product
        end
    end

    % Calculate right-hand side vector b where b(i) = (f, P_i)_{L^2}
    b = zeros(dim, 1);
    for i = 1:dim
        b(i) = sum(f .* T{i});  % L^2 inner product
    end

    % Solve linear system G*coeffs = b for the coefficients
    coeffs = G \ b;

    % Form the projection
    proj = coeffs(1) * T{1};
    for i = 2:dim
        proj = proj + coeffs(i) * T{i};
    end
end