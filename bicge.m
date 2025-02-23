load west0479;
A = west0479;


%Define b so that the true solution is a vector of all ones:
B = full(sum(A,2));



%Set the tolerance and maximum number of iterations:
tol = 1e-6; maxit = 20;


%Use bicg to find a solution at the requested tolerance and number of iterations:
[x0,fl0,rr0,it0,rv0] = bicg(A,B,tol,maxit);

semilogy(0:maxit,rv0/norm(B),'-o');
xlabel('Iteration number');
ylabel('Relative residual');

%%
% using precondition
tic
[L,U] = ilu(A,struct('type','ilutp','droptol',1e-6));
[x1,fl1,rr1,it1,rv1] = bicg(A,B,tol,maxit,L,U);
toc
figure(2)
semilogy(0:it1,rv1/norm(B),'-o');
xlabel('Iteration number');
ylabel('Relative residual');
