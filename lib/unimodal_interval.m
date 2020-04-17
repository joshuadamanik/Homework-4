function [Xa, Xb] = unimodal_interval(X, eps, p, f)
    Xa = X;
    Xb = X;
    p = p / norm(p);

    while true
        lastX = X;
        X = X + eps * p;
        Xb = X;

        eps = 1.5*eps;

        if (f(X) < f(lastX))
            Xa = lastX;
        else
            break;
        end
    end
end