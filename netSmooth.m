function Ft = netSmooth(F, A, alpha, maxeps, n)
%%
%% function Ft = netSmooth(F, A, alpha, maxeps, n)
%%

	
Ft = F;
sa = sum(A);
dsa = diag(1./sa);
An = A*dsa;
counter = 1;
eps = 1+maxeps;
while eps > maxeps && counter < n
    Ft1 = alpha*Ft*An + (1-alpha)*F;
    eps = norm(Ft1-Ft);
    Ft = Ft1;
    counter = counter +1;
end

% Fill the NaNs (for isolated nodes) with original values
Ft(isnan(Ft)) = F(isnan(Ft));

