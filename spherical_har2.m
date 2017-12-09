function [y,y_theta,y_fai]=spherical_har2(m,n,theta,fai)

% this function is for solving the spherical harmonics 
% y is the output 
% n: the order of spherical harmonics
% m: the degree of spherical harmonics
% theta:¸©Ñö½Ç [0,pi]
% fai:·½Î»½Ç [0,2*pi]
% write on 2009.4.7
% edit on 4.12

i=sqrt(-1);
p=legendre(n,cos(theta)); 

% so harmonics function must be cut into [0,m] and [-m,-1].
if m>=0    
    y=sqrt((2*n+1)*factorial(n-m)/4/pi/factorial(n+m))*p((m+1),:).*exp(i*m*fai);
    y_fai=i*m*y;
    y_theta2=-sin(theta)/(1-(cos(theta))^2)*m*cos(theta)*y;

else
    y=(-1)^(-m)*sqrt((2*n+1)*factorial(n+m)/4/pi/factorial(n-m))*p((-m+1),:).*exp(i*m*fai);
    y_fai=i*m*y;
    y_theta2=-sin(theta)/(1-(cos(theta))^2)*(-m)*cos(theta)*y;
    
end
if m>=1
    y_theta1=-sin(theta)/(1-(cos(theta))^2)*sqrt((2*n+1)*factorial(n-m)/4/pi/factorial(n+m))*(n-m+1)*(n+m)*sqrt(1-(cos(theta))^2)*p((m),:).*exp(i*m*fai);
elseif m==0
     if n==0
        y_theta1=0;
     else 
        y_theta1=-sin(theta)/(1-(cos(theta))^2)*sqrt((2*n+1)*factorial(n-m)/4/pi/factorial(n+m))*(n-m+1)*(n+m)*sqrt(1-(cos(theta))^2)*(-1)*factorial(n-1)/factorial(n+1)*p((2),:).*exp(i*m*fai);
     end
else
     y_theta1=-sin(theta)/(1-(cos(theta))^2)*(-1)^(-m)*sqrt((2*n+1)*factorial(n+m)/4/pi/factorial(n-m))*(n+m+1)*(n-m)*sqrt(1-(cos(theta))^2)*p((-m),:).*exp(i*m*fai);
end
   y_theta=y_theta1+y_theta2;


