function [ L ] = laplacian( A )
%LAPLACIAN Calculates the laplacian of the adjacency matrix

D=sum(A,2);
Dd=diag(D);
L=Dd-A;

end

