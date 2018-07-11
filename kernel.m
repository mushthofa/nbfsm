function [result] = kernel(A, method, param1)
%KERNEL calculates graph kernel from the adjacency matrix
%   KCT: commute time kernel
[n,~]=size(A);
result=zeros(n);
[A,OK]=checkData(A);

switch(method)
    case 'CT'
        L=laplacian(A);
        K=pinv(L);
    case 'MD'
        n = size(A,1);
        P=transition(A);
        P(P==0)=1/(n*n);
        t=5; %t=param1;
        I=eye(n);
        Zt=1/t*(I-P)^(-1)*(I-P^t)*P;
        K=Zt*Zt';
    case 'RED'
        n = size(A,1);
        P=transition(A);
        P(P==0)=1/(n*n);
        t=1; %param1;
        I=eye(n);
        Zt=1/t*(I-P)^(-1)*(I-P^t)*P;
        K=Zt*log(Zt')+log(Zt)*Zt';
    case 'RWR'        
        Ds=sum(A,2);
        D=diag(Ds);    
        a=1; %param1;        
        K=inv(D-a*A)*D;
    case 'RW'                
        K=kernel(A,'RWR',0.9999);
    case 'LEXP'
        L=laplacian(A);
        K=expm(-param1*L);
    case 'EXP'
        K=expm(param1*A);
    case 'RL'
        n = size(A,1);
        I=eye(n);
        L=laplacian(A);
        %rad = max(abs(eig(L))); 
        %alf = 1./(2*rad);
        K=inv(I+param1*L); 
    case 'RCT'
        param1=0.9;
        Ds=sum(A,2);
        D=diag(Ds);    
        a=param1;        
        K=pinv(D-a*A);        
    case 'VND'
        n = size(A,1);
        I=eye(n);
        K=inv(I-param1*A);   
    otherwise
        disp('unknown kernel type');
end

result(OK)=K;
result(~OK)=NaN;

    function [A,OK]=checkData(A)
        toKeep=sum(A~=0)>0;
        OK=zeros(size(A));
        A=A(toKeep,toKeep);
        OK(toKeep,toKeep)=1;
        OK=OK==1;
    end
end

