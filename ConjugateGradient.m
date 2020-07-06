%% Conjugate Gradient Method
clear; clc;
fprintf('Performing Conjugate Gradient Method...\n')
syms x y a b
% example in class from Lecture 5 June 9 - page 4
f_handle = @(x,y) 5*x.^2 + y.^2 + 4*x*y -14*x - 6*y + 20;
f = 5*(x.^2) + (y.^2) + (4*x*y) -(14*x) - (6*y) + 20;
x0 = 0;
y0 = 10;
eps = 1e-6;

f_handle = @(x,y) 100*(y - x).^2 + (1 - x).^2;
f = 100*(y-x).^2 + (1-x).^2;
x0 = -1.2;
y0 = 1;
eps = 0.00005;

initial_guess = [ x0 ; y0 ];

%Step 2 - Calculate the Negative Gradient - towards the minimum
delta_f(x,y) = - gradient(f);
s0 = double(delta_f(x0,y0));
f_handle(initial_guess(1),initial_guess(2));

% Try plotting
figure(1); clf;
subplot 111    % Two rows One Column Second Place
[x,y] = meshgrid(-50:1:50);
f1 = 100*(y-x).^2 + (1-x).^2;
surf(x,y,f1);
hold on

iter = 1;
condition = eps+1;

fid = fopen('ConjugateGradient_vals.csv','wt');
fprintf(fid,'iteration,sx,sy,x,y,f(x y),f_change\n');

while(condition > eps)
    fprintf('---------\nIteration %d\n',iter)
    f_prev = f_handle;
    %Step 3 - Search Direction #1 Iteration
    if(iter == 1)
        x1 = initial_guess + a*s0;
        f_handle_x1 = f_handle(x1(1),x1(2));
        d = diff(f_handle_x1, a);
        b = vpa(double(solve(d, a)));
        x = vpa(double(initial_guess + b.*s0));
        f_eval = f_handle(x(1),x(2));
        plot(x(1),x(2),'ro','LineWidth',5)
    %Step 4 - Search Direction #N Iteration
    else
        %Step 5 - Go back to Step 3 and calculate the values of x1
        %   replace a with b
        delta_f2 = vpa(double(delta_f(x(1),x(2))));
        reshape(delta_f2,1,2);
        beta = (((delta_f2(1)).^2) + (delta_f2(2)).^2 )/...
            ((s0(1)).^2 + (s0(2)).^2);
        s1 = delta_f2 + beta.*s0;
        xn = x + a.*s1;
        f_handle_xn = f_handle(xn(1),xn(2));
        d1 = diff(f_handle_xn, a);
        b2 = vpa(double(solve(d1, a)));
        xn = x + b2.*s1;
        f_eval = f_handle(xn(1),xn(2));
        % dif b/w prev & curr
        condition = abs(double(f_eval - f_prev_eval));
        plot(xn(1),xn(2),'ro','LineWidth',5)
        x = xn;
        s0 = s1;
    end
    % print to file
    fprintf(fid,'%d,%f,%f,%f,%f,%f,%e\n',...
        iter,s0(1),s0(2),x(1),x(2),f_handle(x(1),x(2)),condition);
    
    % update vals
    f_prev_eval = f_eval;
    iter = iter + 1;
    
    fprintf('s direction: [%f, %f]\t',s0(1),s0(2))
    fprintf('x = [%f, %f]\t\t f(x) = %e\n', ...
        x(1),x(2), f_handle(x(1),x(2)))
    fprintf('f_change = %e\n',condition)
    
end
plot(xn(1),xn(2),'go','LineWidth',5)
fprintf('---------\nAfter %d iterations, ',iter-1)
fprintf('the Minimizer is located at f(x1 = %f, x2 = %f) = %f\n',...
        xn(1),xn(2),f_handle(xn(1),xn(2)));
fclose(fid);
