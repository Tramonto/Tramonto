function [x2,flag,relres,iter,resvec] = run_gmres(A,x,b)
tol = 1e-3;  maxit = 1; restart=200;
%[L,U,P]=luinc(A,'0'); %luinc
%D=diag(A); % Jacobi scaling
[x2,flag,relres,iter,resvec] = gmres(@afun,b,restart,tol,maxit,@mfun,[],x);
x1 = x2;
    function y = afun(x)
        y = A*x;
    end
 
%    function y = mfun(r) %luinc
%        y=L\(P*r);
%        y=U\y;
%    end
%    function y = mfun(r) %Jacobi scaling
%        y=r./D;
%    end
    function y = mfun(r) %No Prec
        y=r;
    end
end