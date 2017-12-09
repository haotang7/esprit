function [x]=gauss(a,b)
n=length(a);
x=zeros(n,1);
a=[a b];
for k=1:n-1
    max=k;
    for i=k+1:n
        if a(i,k)>a(max,k)
            max=i;
        end
    end
    temp=a(k,k:n+1);
    a(k,k:n+1)=a(max,k:n+1);
    a(max,k:n+1)=temp;
    for i=k+1:n
        a(i,k)=-a(i,k)/a(k,k);
        a(i,k+1:n+1)=a(i,k+1:n+1)+a(i,k)*a(k,k+1:n+1);
    end
end
x(n,1)=a(n,n+1)/a(n,n);
for i=n-1:-1:1
    sum=0;
    for j=i+1:n
        sum=sum+x(j,1)*a(i,j);
    end
    x(i,1)=(a(i,n+1)-sum)/a(i,i);
end