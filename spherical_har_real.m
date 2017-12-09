function y=spherical_har_real(m,n,theta,fai)

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
% Because legendre function is just for [0,m], 
% so harmonics function must be cut into [0,m] and [-m,-1].
if m>0    
    y=sqrt((2*n+1)*factorial(n-m)/4/pi/factorial(n+m))*p((m+1),:).*exp(i*m*fai);
    y=sqrt(2)*real(y);
else if m<0
        m=-m;
        y=sqrt((2*n+1)*factorial(n-m)/4/pi/factorial(n+m))*p((m+1),:).*exp(i*m*fai);
        y=sqrt(2)*imag(y);
    else
         y=sqrt((2*n+1)/4/pi)*p(1,:);
    end
end