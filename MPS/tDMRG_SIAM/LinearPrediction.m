function [x_LP] = LinearPrediction(x,Nmax)
%Linear prediction to extrapolates Chebyshev moments / real-time correlations or currents 

%Parameters for extrapolation:
w=1; %weigthing function (usually chosen to be constant)
delta = 1e-10; %regularization parameter for pseudoinverse (usually ranging between 1E-4 - 1E-10)
Nc=length(x); %number of hardcore computed moments / time points
nfit=floor(1*Nc/4); %defines interval of data points used for fit of cost function: {Nc-nfit,...,Nc}; used to eliminate short-time behavior
p=min(floor((Nc-nfit)/2),120); %number of points used in extrapolations 

%IMPORTANT: check robustness of output by changing extrapolation
%parameters! If small parameter changes lead to drastic deviations in the extrapolated
%data, linear prediction might NOT be applicable!

%safeguard
if size(x,2) ~=1 && size(x,1) ~=1
   fprintf('\n ERR: input not a vector, dim(x) = %g \n', size(x));
end

if size(x,2) ~= 1 && size(x,1) == 1
   x=x.';
end
   
%set up normal equations to minimize cost function: R*a = -r
R=zeros(p,p);
r=zeros(p,1);
for k=1:p
    for l=1:p
        for m=(Nc-nfit+1):Nc
            R(k,l) = R(k,l)+w*conj(x(m-k))*x(m-l);
        end
    end
    for m=(Nc-nfit+1):Nc
        r(k)= r(k)+w*conj(x(m-k))*x(m);
    end
end

%solve normal equation using pseudo inverse: a = -R^-1 r
a =-pinv(R,delta)*r;

%form companion matrix
M=diag(ones(p-1,1),-1);
for k=1:p
    M(1,k) = - a(k);
end

%eleminate eigenvalues > 1
[U,Lam] = eig(M);
for k=1:p
    if abs(Lam(k,k))>1
        Lam(k,k)=0;
    end
end
UI=inv(U);

%generate 'truncated' matrix M
M2 = U*Lam*UI;


%generate predicted moments/correlators
temp = (x(Nc:-1:Nc-p+1));
x_LP=x; %first Nc moments in x_LP are the computed ones
for k=Nc+1:Nmax
    temp = M2*temp;
    x_LP(k)=temp(1);
end


end%%function
