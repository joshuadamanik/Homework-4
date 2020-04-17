function g = grad_central_diff(X, eps, f)
    dim = length(X);
    E = eye(dim);
    
    g = zeros(dim, 1);
    for i = 1:dim
        g(i) = (f(X + eps*E(:,i)) - f(X - eps*E(:,i))) / (2*eps);
    end
end