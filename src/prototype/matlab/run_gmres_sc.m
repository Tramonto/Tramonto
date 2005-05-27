function [xg,flag,relres,iter,resvec] = run_gmres_sc(A,x,b)
tol = 1e-3;  maxit = 1; restart=200;
nbeads=18;
nunk=40;
n=161;
nn=nunk*n;
n1=36*n; % n1=size of 1,1 block
[prow,pcol]=makenewpnew(); 
Ap=A(prow,pcol);
xp=x(pcol);
bp=b(prow);
A11=Ap(1:n1,1:n1);
A22=Ap(n1+1:nn,n1+1:nn);
A21=Ap(n1+1:nn,1:n1);
A12=Ap(1:n1,n1+1:nn);
x1=xp(1:n1);
x2=xp(n1+1:nn);
b1=bp(1:n1);
b2=bp(n1+1:nn);
bs=b2-A21*A11invfun(b1);
%D=diag(A22); %Jacobi scaling
%[L,U,P]=lu(A22); %lu
%[L,U,P]=luinc(A22,'0'); %luinc
[x2,flag,relres,iter,resvec] = gmres(@afun,bs,restart,tol,maxit,@mfun);
xp(n1+1:nn) = x2;
btmp=b1-A12*x2;
xp(1:n1)= A11invfun(btmp);
xg=xp(pcol);
    function y2 = afun(x2)
        xtmp=A12*x2;
        xtmp=A11invfun(xtmp);
        y2 = A22*x2 - A21*xtmp;
    end
 
    function y2 = mfun(r2) %no prec
        y2=r2;
    end
 %   function y2 = mfun(r2) %Jacobi scaling
 %       y2=r2./D;
 %   end
 %   function y2 = mfun(r2) %luinc or lu
 %       y2=L\(P*r2);
 %       y2=U\y2;
 %   end
    function z = A11invfun(w)
        z=w; % init z
        z(1:n)=A11(1:n,1:n)\z(1:n); % Do first solve separately (no matvec needed)
        for i=2:2*nbeads
            nstrt=(i-1)*n+1;
            nstop=nstrt+n-1;
            xtmp=z(1:nstrt-1); % get the known solution values
            Amv=A11(nstrt:nstop,1:nstrt-1);
            btmp=z(nstrt:nstop) - Amv*xtmp;
            z(nstrt:nstop)=A11(nstrt:nstop,nstrt:nstop)\btmp;
        end
    end
end