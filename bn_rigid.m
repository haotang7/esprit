function bn=bn_rigid(kr,ka,n)

% This is for calculating bn of rigid sphere;
% bn is the output;
% k: wave number for frequency w and the speed of sound c, k=w/c
% r: the radius of definition
% a: the radius of a rigid sphere
% n: the order of spherical harmonics fuction
% written on 2009.8.23

i=sqrt(-1);
j_kr=sqrt(pi/2./kr).*besselj(0.5+n,kr);     %the first spherical bessel function
y_kr=sqrt(pi/2./kr).*bessely(0.5+n,kr);     %the second spherical bessel function
% the first hankel function
h_kr=j_kr+i*y_kr;
% j_ka=sqrt(pi/2./ka).*besselj(0.5+n,ka);     %the first spherical bessel function
% y_ka=sqrt(pi/2./ka).*bessely(0.5+n,ka);     %the second spherical bessel function
% from the reference P550, 
% get dj/dkr=(n*j(n-1,kr)-(n+1)*j(n+1,kr))/(2*n+1);  
% take j1=dj/dkr
j1_ka=(n*sqrt(pi/2./ka).*besselj(-0.5+n,ka)-(n+1)*sqrt(pi/2./ka).*besselj(1.5+n,ka))/(2*n+1);
% take y1=dy/dkr
y1_ka=(n*sqrt(pi/2./ka).*bessely(-0.5+n,ka)-(n+1)*sqrt(pi/2./ka).*bessely(1.5+n,ka))/(2*n+1);
% the first hankel function
% h_ka=j_ka+i*y_ka;
% take h1=dh/dkr
h1_ka=j1_ka+i*y1_ka;
% bn of rigid sphere is
bn=4*pi*i^n*(j_kr-j1_ka./h1_ka.*h_kr);