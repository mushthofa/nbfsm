function kernmat = kernelNorm(kernmat)
% Normalise kernel matrix
dd = diag(kernmat);
ddmat = dd*dd';
ddmat = sqrt(ddmat);
kernmat = kernmat./ddmat;

end