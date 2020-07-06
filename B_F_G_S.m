%% BFGS Method 
clear; clc;
fprintf('Performing BFGS Method...\n')
% Exercise
% f_handle = @(x,y) 12*x.^2 + 4*y.^2 - 12*x*y + 2*x;
% f =  12*x.^2 + 4*y.^2 - 12*x*y + 2*x;
% x0 = -1;
% y0 = -2;
% eps = 0.00005;

initial_guess=[1;2];
xprev = initial_guess;
eps = 0.00005;
condition = eps + 1;
iter = 1;
Bprev = eye(2);   % B-1 = identity matrix 2x2

f_handle=@(x,y) x.^4 - 2.*y.*x.^2 + y.^2 + x.^2 - 2.*x + 5;

% Actual functions
syms x y
f = x.^4 - 2*y.*x.^2 + y.^2 + x.^2 - 2*x + 5;
f_delta = gradient(f);

% Set alpha to ad hoc value; we're not optimizing it here
alpha = 1.0;
fid = fopen('BFGS_vals.csv','wt');
fprintf(fid,'iteration,sx,sy,x,y,f(x y),f_change\n');

while condition > eps
    fprintf('---------\nIteration %d\n',iter)
    s = -(Bprev)*(subs(f_delta,[x y],[xprev(1) xprev(2)]));
    xk = vpa(xprev + alpha*s);
    condition = abs(f_handle(xk(1),xk(2)) - ...
                    f_handle(xprev(1),xprev(2)));
    d = xk - xprev;
    yn = subs(f_delta,[x y],[xk(1) xk(2)]) - ...
            subs(f_delta,[x y],[xprev(1) xprev(2)]);
    
    % Update using Sherman-Morrison Formula
    Bk = Bprev-(Bprev*yn*d.'+d*yn.'*Bprev)/(d.'*yn) + ...
            (1+(yn.'*Bprev*yn)/(d.'*yn))*(d*d.')/(d.'*yn);
    
    % print to file
    fprintf(fid,'%d,%f,%f,%f,%f,%f,%e\n',...
        iter,s(1),s(2),xk(1),xk(2),f_handle(xk(1),xk(2)),condition);
    
    xprev = xk;         % update xprev before proceeding
    Bprev = Bk;         % update Bprev before proceeding
    iter = iter + 1;
    
    fprintf('s direction: [%f, %f]\t',s(1),s(2))
    fprintf('x = [%f, %f]\t\t f(x) = %f\n', ...
        xk(1),xk(2), f_handle(xk(1),xk(2)))
    fprintf('f_change = %e\n',condition)
%     csvwrite('BFGS_vals.csv',s(1),s(2),xk(1),xk(2),f_handle(xk(1),xk(2)))
end
fprintf('---------\nAfter %d iterations, ',iter-1)
fprintf('the Minimizer is located at f(x1 = %f, x2 = %f) = %f\n',...
        xk(1),xk(2),f_handle(xk(1),xk(2)));
fclose(fid);
