# Minimization-Optimization-Methods
Author: Mohammad Einaam Alim, Raphael Barata de Oliveira
Various methods for optimization

## Description
The problem is finding the minimum of the function f(x)=x^4+2xy^2+y^2-2x+5, using the Broyden-Fletcher-Goldfarb-Shanno (BFGS) Method with an initial value (x_0,y_0)=(1,2)^T. The iteration stops when the estimated error value |f(x_(n+1))-f(x_n)|<0.00005. The same tolerance is used for the remaining 4 subproblems in #2 and #3 which each respectively use Conjugate Gradient Method, Steepest Descent Method, SQP Method, and Reduced Gradient Method.

## Discussion
Each program initializes the initial guess and function/function-handle, epsilon tolerance, and initial variables. Each uses a loop with a break-out condition for when the estimated error value is less than the specified epsilon tolerance parameter, which shows the precise minimizer is found. <br>
The program implements BFGS Method by first initializing (x_0,y_0)=(1,2)^T. Then it proceeds to a loop that iteratively obtains the gradient of f(x_0) and the search direction s_n=B_k ∇f(x_n) and then determines x_(n+1)=x_n+αx_(n+1) s_n. The Sherman-Morrison Formula is used to determine each iteration’s respective B_k and stores previous values for the next iteration’s use. Alpha is set to an ad hoc value of 1.0, as we are not optimizing it here. <br>

For 2A, the program implements Conjugate Gradient Method by first initializing (x_0,y_0)=(-1.2,1.0)^T and f(x)=100(y-x)^2+(1-x)^2. The program determines the first search direction  s_0=-∇f(x_0,y_0) and optimizes alpha. For following iterations, β was obtained and used to determine the search direction s_(n+1)=-∇f(x_n,y_n)+βs_n and x_(n+1)=x_n+αs_n. <br>

For 2B, the program implements Steepest Descent Method by first initializing (x_0,y_0)=(-1.2,1.0)^T and f(x)=100(y-x)^2+(1-x)^2. Similarly, it shares the same first iteration step as Conjugate Gradient Method with the search direction being s_0=-∇f(x_0,y_0); however, for all of the following iterations, the method updates the search direction with the following relationship: s_n=-∇f(x_n,y_n) and x_(n+1)=x_n+αs_n. <br>

For 3A, the program implements SQP Method by first initializing (x_0,y_0)=(0.0,√2)^T and f(x)=100(y-x)^2+(1-x)^2 subject to x^2+y^2=2. For each iteration, a matrix A and B are calculated using the Lagrangian and updating a value λ by an obtained value of v found by taking the third element of P. The matrix P is solved for by using the matrices A and B:
[■(〖∇^2〗_xx L&-∇h@-(∇h)&0)][■(P_k@Y_k )]=[■(〖-∇〗_x L@h)]
in the form of AP=B. <br>

For 3B, the program implements Reduced Gradient Method by first initializing (x_0,y_0)=(≈0.0,√2)^T and f(x)=100(y-x)^2+(1-x)^2 subject to x^2+y^2=2. This approach uses the so-called ‘hard’ way that uses explicit elimination. Similar steps to SQP Method were taken, except this program instead sets and finds a value γ (r in the code) that determines the search direction P_k=Z_k P_z+Y_k P_γ and updates the iteration value x_(k+1)=x_k+P_k.
