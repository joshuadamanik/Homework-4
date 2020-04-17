function X_star = fibonacci_search(Xa, Xb, k0, f)
    fib0 = Fib(k0);
    
    k = k0 - 2;
    fib_min = 0;
    fib_max = fib0;
    
    while k > 0
        del_fib = Fib(k);
        fib_range = [fib_min, fib_min + del_fib, fib_max - del_fib, fib_max];
        
        X_range = Xa + (Xb-Xa) * fib_range / fib0;
        for i = 1:4
            Z_range(i) = f(X_range(:,i));
        end
        [~, imin] = min(Z_range);
        
        X_star = X_range(:,imin);
        
        if (imin < 3)
            fib_max = fib_range(3);
        else
            fib_min = fib_range(2);
        end
        k = k - 1;
    end

    X = X_star;
end