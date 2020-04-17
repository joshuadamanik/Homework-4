function [Q, R, S, Ak, Bk] = cost_function_param(N, a, c)
    dt = 20 / N;
    
    y0 = 1000;
    vy0 = 0;
    
    Y0 = [y0, vy0]';
    A = [1, dt;0, 1];
    B = [0, dt]';
    
    Ak(1:N,:) = {A};
    Bk_col(1:N,1) = {B};
    Bk(1:N,1:N) = {zeros(2,1)};
    
    for k = 2:N
        for j = k:N
            Ak{j,:} = Ak{j,:} * A;
            Bk_col{j,:} = A*Bk_col{j,:};
        end
    end
    
    for k = 1:N
        for j = k:N
            Bk(j,k) = Bk_col(j-k+1);
        end
    end
    
    Ak = cell2mat(Ak);
    Bk = cell2mat(Bk);
    
    P = [];
    for i = 1:N
        Pi = [a, 0;0, 0];
        if i == N
            Pi = Pi + [c, 0;0, 0];
        end
        P = blkdiag(P,Pi);
    end

    Q = Ak'*P*Ak;
    R = Bk'*P*Bk + eye(N);
    S = Ak'*P*Bk;
end