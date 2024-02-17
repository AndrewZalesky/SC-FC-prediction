function xmat=map_to_matrix(x)

N=size(x,1); %number of subjects
J=size(x,2); %number of connections

%Reconstruct matrices
K=(2*J + 1/4)^(1/2) + 1/2; %number of nodes
ind_upper=find(triu(ones(K,K),1));
xmat=zeros(K,K,N); ymat=zeros(K,K,N);
for i=1:N
    tmp=zeros(K,K);
    tmp(ind_upper)=x(i,:);
    tmp=tmp+tmp';
    tmp=tmp+eye(K);
    %xmat(:,:,i)=logm(tmp);
    xmat(:,:,i)=tmp;
end

