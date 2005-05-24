function [x2,flag,relres,iter,resvec] = run_gmres(A,x,b)
tol = 1e-12;  maxit = 1; restart=200;
[L,U,P]=luinc(A,'0');
[x2,flag,relres,iter,resvec] = gmres(@afun,b,restart,tol,maxit,@mfun,[],x);
x1 = x2;
    function y = afun(x)
        y = A*x;
    end
 
    function y = mfun(r)
        y=L\(P*r);
        y=U\y;
    end
end