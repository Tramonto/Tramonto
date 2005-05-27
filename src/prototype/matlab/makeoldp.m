n=161;
nunk=40;
pold1 = [1:161:6440];
pold = zeros(40*n,1);
for i=1:161
    pold((i-1)*nunk+1:i*nunk) = pold1+(i-1);
end