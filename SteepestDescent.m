%% Steepest Descent
clear; clc;
fprintf('Performing Steepest Descent Method...\n')
syms x y a b

% example in class from Lecture 5 June 9 - page 4
f_handle = @(x,y) 5*x.^2 + y.^2 + 4*x*y -14*x - 6*y + 20;
f = 5*(x.^2) + (y.^2) + (4*x*y) -(14*x) - (6*y) + 20;
x0 = 0;
y0 = 10;
eps = 1e-6;

f_handle = @(x,y) 100*(y - x).^2 + (1 - x)^2;  % TYPO HERE
f = 100*(y-x).^2 + (1-x)^2;
x0 = -1.2;
y0 = 1;
eps = 0.00005;

% Confirming SQP results for old problem
% f = (x.^2+y.^2);
% x0 = 0;
% y0 = sqrt(2);
% % c = (x.^2)+(y.^2) - 2;              % constraints h(xk)
% f_vector = @(x,y) (x^2+y^2);
initial_guess = [ x0 ; y0 ];

%Step 2 - Calculate the Negative Gradient - towards the minimum
delta_f(x,y) = - gradient(f);
s0 = double(delta_f(x0,y0));
f_handle(initial_guess(1),initial_guess(2));

iter = 1;
condition = eps+1;

% Try plotting
figure(1); clf;
subplot 111    % Two rows One Column Second Place
[x,y] = meshgrid(-50:1:50);
f1 = 100*(y-x).^2 + (1-x.^2);
% f1 = (x.^2+y.^2);
surf(x,y,f1);
hold on

fid = fopen('SteepestDescent.csv','wt');
fprintf(fid,'iteration,sx,sy,x,y,f(x y),f_change\n');

while(condition > eps)
    fprintf('---------\nIteration %d\n',iter)
    f_prev = f_handle;
    %Step 3 - Search Direction #1 Iteration
    if(iter == 1)
        x1 = initial_guess + a*s0;
        f_vector_x1 = f_handle(x1(1),x1(2));
        d = diff(f_vector_x1, a);
        b = vpa(double(solve(d, a)));
        xn = vpa(double(initial_guess + b.*s0));
        f_eval = f_handle(xn(1),xn(2));
        s = s0;
        plot(xn(1),xn(2),'ro','LineWidth',5)
    %Step 4 - Search Direction #n Iteration
    else
        s1 = vpa(double(delta_f(xn(1),xn(2))));
        x2 = xn + a.*s1;
        f_vector_x2 = f_handle(x2(1),x2(2));
        d = diff(f_vector_x2, a);
        b = vpa(double(solve(d, a)));
        xn = vpa(double(xn + b.*s1));
        f_eval = f_handle(xn(1),xn(2));
        condition = abs(double(f_eval - f_prev_eval));
        s = s1;
        plot(xn(1),xn(2),'ro','LineWidth',5)
    end
    
    % print to file
    fprintf(fid,'%d,%f,%f,%f,%f,%f,%e\n',...
        iter,s(1),s(2),xn(1),xn(2),f_handle(xn(1),xn(2)),condition);
    
    % values to hold previous values
    f_prev_eval = f_eval;
    iter = iter + 1;
    
    fprintf('s direction: [%f, %f]\t',s(1),s(2))
    fprintf('x = [%f, %f]\t\t f(x) = %e\n', ...
        xn(1),xn(2), f_handle(xn(1),xn(2)))
    fprintf('f_change = %e\n',condition)
end

fprintf('---------\nAfter %d iterations, ',iter-1)
fprintf('the Minimizer is located at f(x1 = %f, x2 = %f) = %f\n',...
        xn(1),xn(2),f_handle(xn(1),xn(2)));
fclose(fid);
