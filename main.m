%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOMEWORK #4
% Joshua Julian Damanik (20194701)
% AE551 - Introduction to Optimal Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all;
addpath('lib');

X0 = [1000, 0]';
N = 50;

U = zeros(N,1);
eps = 0.01;
k = 15;

a_list      = [   1,   10,  0.1,    0,     0,    0];
c_list      = [   1,    1,    1,    1, 0.005, 0.01];
marker_list = {'o-', '+-', '^-', 'v-',  'x-', '.-'};
legends = {};

plot_data = zeros(N+1,length(a_list),3);

fprintf('U* program\tU* analytic\tError\n');
for j = 1:length(a_list)
    a = a_list(j);
    c = c_list(j);
    
    [Q, R, S, Ak, Bk] = cost_function_param(N, a, c);
    f = @(U) 0.5*(X0'*Q*X0 + U'*R*U + 2*X0'*S*U);
    %f = @cost_function;
    
    quasi = quasi_newton_class(eye(N));
    
    for n = 1:100
        p = quasi.bgfs(U, eps, f);
        [Xa, Xb] = unimodal_interval(U, eps, p, f);
        U_star = fibonacci_search(Xa, Xb, k, f);
        
        err = norm(U_star - U);
        if (err < eps * 0.01)
            break;
        end
        U = U_star;
    end
    
    U_anal = -R\S'*X0;
    Erms = norm(U_anal-U_star)/norm(U_anal)*100;
    fprintf('%.4f\t%.4f\t%.8f\n', norm(U_star), norm(U_anal), Erms);
    
    X = Ak*X0+Bk*U_star;
    
    plot_data(1,j,1) = X0(1);
    plot_data(1,j,2) = X0(2);
    for n = 1:N
        plot_data(n+1,j,1) = X(2*n-1);
        plot_data(n+1,j,2) = X(2*n);
    end
    
    plot_data(:,j,3) = [U_star; 0];
    
    legends = [legends, {sprintf('a=%.1f, c=%.3f', a, c)}];
end

title = {'Vertical Position', 'Vertical Speed', 'Vertical Input'};
ylabels = {'Y (m)', 'Vy (m/s)', 'U (m/s2)'};
legend_position = {'NorthEast', 'SouthEast', 'SouthEast'};
for i = 1:3
    figure('Name', title{i}), hold on;
    set(gcf, 'Position', [100*i, 100*i, 1000, 400]);
    for j = 1:length(a_list)
        figure(i),plot(5000 / N * (0:N), plot_data(:,j,i), marker_list{j});
    end
    grid on;
    xlabel('x (m)');
    ylabel(ylabels{i});
    legend(legends, 'Location', legend_position{i});
end

%% Functions

function [J, Y, Vy] = cost_function(U)

N = size(U,1);
dt = 20 / N;

a = 1;
c = 1;

Y = zeros(N, 1);
Vy = zeros(N,1);

Y(1) = 1000;
Vy(1) = 0;

for i = 2:(N+1)
    Vy(i) = Vy(i-1) + U(i-1)*dt;
    Y(i) = Y(i-1) + Vy(i-1)*dt;
end

J = c./2.*Y(N+1).^2 + 1./2.*sum(a.*(Y(2:N+1).^2) + U(1:N).^2);
end

