function [prow,pcol]= makenewpnew()

n=161;
p1 = [1:n];
p2 = [1:10];
p2c = [11:n n];
p3 = [1:12];
p3c = [13:n];
prow = zeros(40*n,1);
pcol = prow;
prow = [ 4*n+p1 6*n+p1 8*n+p1 10*n+p1 12*n+p1 14*n+p1 16*n+p1 18*n+p1 20*n+p1 22*n+p1 24*n+p1 26*n+p1 28*n+p1 30*n+p1 32*n+p1 34*n+p1 36*n+p1 38*n+p1 39*n+p1 37*n+p1 35*n+p1 33*n+p1 31*n+p1 29*n+p1 27*n+p1 25*n+p1 23*n+p1 21*n+p1 19*n+p1 17*n+p1 15*n+p1 13*n+p1 11*n+p1 9*n+p1 7*n+p1 5*n+p1   p2 n+p3   2*n+p2c 3*n+p3c   2*n+p2 3*n+p3   p2c     n+p3c     ];
pcol = [ 4*n+p1 6*n+p1 8*n+p1 10*n+p1 12*n+p1 14*n+p1 16*n+p1 18*n+p1 20*n+p1 22*n+p1 24*n+p1 26*n+p1 28*n+p1 30*n+p1 32*n+p1 34*n+p1 36*n+p1 38*n+p1 39*n+p1 37*n+p1 35*n+p1 33*n+p1 31*n+p1 29*n+p1 27*n+p1 25*n+p1 23*n+p1 21*n+p1 19*n+p1 17*n+p1 15*n+p1 13*n+p1 11*n+p1 9*n+p1 7*n+p1 5*n+p1   p2 n+p3   p2c     n+p3c     2*n+p2 3*n+p3   2*n+p2c 3*n+p3c   ];
n1=36*n;
end