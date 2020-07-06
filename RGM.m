%% RGM Method
% Implementation that utilizes EXPLICIT Elimination
clear; clc;
fprintf('Performing Reduced-Gradient Method...\n')
syms x y a b ri

% Example problem
f = exp(3*x+4*y);
x0 = -0.7;
y0 = -0.7;
c = (x.^2)+(y.^2) - 2;              % constraints h(xk)
f_handle = @(x,y) exp(3*x+4*y);

% Homework problem
% f = (x.^2+y.^2);                    % original
f = 100*((y-x).^2)+((1-x).^2);
f_handle = @(x,y) 100*((y-x).^2)+((1-x).^2);
% Tweaking the initial guess was said to be allowed, and for this specific
%   system, x0 = 0 yields catastrophic results for the general solution.
%   Therefore, we will consider x0 to be very small and close to zero.
x0 = 1e-6;
y0 = sqrt(2);
c = (x.^2)+(y.^2) - 2;              % constraints h(xk)
% % f_handle = @(x,y) (x^2+y^2);      % original
% f_handle = @(x,y) 100*((y-x).^2)+((1-x).^2);

% Try plotting
% figure(1); clf;
% subplot 111    % Two rows One Column Second Place
% [s,t] = meshgrid(-3:0.5:3);
% f1 = 100*((t-s).^2)+((1-s).^2);
% surfc(s,t,f1,'FaceAlpha',0.5);
% hold on
% grid on
% f1 = s.^2+t.^2 - 2;
% surfc(s,t,f1,'FaceAlpha',0.25);


% Initial Guess values
eps = 0.00005;
lambda = -0.018616;
xprev = [x0; y0];       % initial guess
r = 0;                  % call gamma r bc it looks like it

% Initial Step - Gradient f
delta_f = gradient(f);
delta_2_f_x = gradient(delta_f(1));
delta_2_f_y = gradient(delta_f(2));
% delta_2_f = [delta_2_f_x delta_2_f_y]

% Use identity matrix to get second gradient f
%   We need a conditional for when we don't have y in the 2nd gradient,
%    then we need to account for those components in the gradient f

i = 0;
if((numel(delta_2_f_x)==1) && (numel(delta_2_f_y)==1))
    A = eye(2);
    delta_2_f = [delta_2_f_x delta_2_f_y] .* A;
else
    % normal condition, where we have y components present
    % pump on, valve off
    delta_2_f = [delta_2_f_x delta_2_f_y];
    i = 1;
end

% Initial Step - Gradient h
delta_h = gradient(c);
delta_2_h_x = gradient(delta_h(1));
delta_2_h_y = gradient(delta_h(2));

% Use identity matrix to get second gradient h
%   We need a conditional for when we don't have y in the 2nd gradient,
%    then we need to account for those components in the gradient h

j = 0;
if((numel(delta_2_h_x)==1) && (numel(delta_2_h_y)==1))
    A = eye(2);
    delta_2_h = [delta_2_h_x delta_2_h_y] .* A;
else
    % normal condition, where we have y components present
    % pump on, valve off
    delta_2_h = [delta_2_h_x delta_2_h_y];
    j = 1;
end

iter = 1;
xk = xprev;
condition = 1;          % make a do-while loop

fid = fopen('RGM_vals.csv','wt');
fprintf(fid,'iteration,sx,sy,x,y,f(x y),f_change\n');

while (condition > eps)
    fprintf('---------\nIteration %d\n',iter)
    
    % Calculate gradient f
    delta_f_0 = vpa(double(subs(delta_f,[x y],[xprev(1) xprev(2)])));
    delta_2_f_0 = (double(subs(delta_2_f,[x y],[xprev(1) xprev(2)])));
    delta_h_0 = vpa(double(subs(delta_h,[x y],[xprev(1) xprev(2)])));
    
    % if we have symbolic variables present, plug in values
    %   (get rid of them)
    if(i == 1)
        delta_2_f = vpa(double(subs(delta_2_f,[x y],[xprev(1) xprev(2)])));
    end
    if(j == 1)
        delta_2_h = vpa(double(subs(delta_2_h,[x y],[xprev(1) xprev(2)])));
    end
    
    % Obtain gradient xL
    delta_x_L = delta_f_0 - lambda*delta_h_0;
    delta_2_xx_L = vpa(delta_2_f_0 - lambda*delta_2_h);
    
    % this is where reduced method varies from SQP
    delta_t_h = delta_h.';
    z_k = [(-inv(delta_t_h(1))*(delta_t_h(2))); 1];
    y_k = (1/(sqrt(x.^2 + y.^2)))*[x;y];
    z_k_0 = vpa(double(subs(z_k,[x y], [xprev(1) xprev(2)])));
    y_k_0 = vpa(double(subs(y_k,[x y], [xprev(1) xprev(2)])));
    
    p_z = (-(z_k_0).'*(delta_x_L))/((z_k_0).'*delta_2_xx_L*z_k_0);
    
%     x_k = xk + p_z - [y_k_0(1); y_k_0(2)]*ri
    x_k = xk + z_k*p_z;
    
    
    phi = ((x_k(1) + ri*(y_k(1))).^2) + ((x_k(2) + ri*(y_k(2))).^2) - 2;
    phi_diff = diff(phi,ri);
    
    rel = phi/phi_diff;
    phi_new = vpa(subs(phi,[x y ri],[xk(1) xk(2) r]));
    
    % find next r
    rnew = ri - ((x_k(1) + ri*y_k(1)).^2 + (x_k(2) + ri*y_k(2)).^2 -2)/...
        (2* (y_k(1)*((x_k(1)) + ri*y_k(1))) + ...
        (2* (y_k(2)*((x_k(2)) + ri*y_k(2)))));
    rnew = vpa(subs(rnew,[x y ri],[xk(1),xk(2) r]));
    r = rnew;
    lambda = lambda + r;
    % Search direction, p0
    p0 = vpa(z_k_0*p_z) + (y_k_0*rnew);
    xk = xk + p0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('RNEW = %f\n',rnew)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    condition = abs(f_handle(xk(1),xk(2))-f_handle(xprev(1),xprev(2)));
    
    % print to file
    fprintf(fid,'%d,%f,%f,%f,%f,%f,%e\n',...
        iter,p0(1),p0(2),xk(1),xk(2),f_handle(xk(1),xk(2)),condition);
    
    % update xprev before proceeding
    xprev = xk;
    iter = iter + 1;

    fprintf('P direction: [%f, %f]\t',p0(1),p0(2))
    fprintf('x = [%f, %f]\t\t f(x) = %e\n', ...
        xk(1),xk(2), f_handle(xk(1),xk(2)))
    fprintf('f_change = %e\n',condition)
end

fprintf('---------\nAfter %d iterations, ',iter-1)
fprintf('the Minimizer is located at f(x1 = %f, x2 = %f) = %f\n',...
        xk(1),xk(2),f_handle(xk(1),xk(2)));
fclose(fid);
