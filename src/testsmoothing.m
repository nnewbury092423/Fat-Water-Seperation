clear
close all
clc
lambda = .24; 
theta = .1;
stepsize =1;
Mask = ones(10,10);
phi = rand(10,10);
W = ones(10,10);
F = ones(10,10);
val = (W+F*exp(1i*theta));
S_exp = abs(val)*exp(.1i);
S_exp(5,5) = abs(val(5,5))*exp(.5i);
for n =1:10000
   dphi = fatwatfuncs.calcdphi(W, F, theta, phi, S_exp,lambda, zeros(size(S_exp)),Mask);
   %dphi = fatwatfuncs.calcdphi(W, F, theta, phi, S_exp,lambda, dphi,Mask);
   phi = phi + (stepsize)*dphi;
   %funcs = fatwatfuncs.Psi(W,F,theta,phi,S_exp,lambda,Mask);
   xgrad = gradient(phi);
   %[xxgrad, ~] = gradient(xgrad);
   %[~, yygrad] = gradient(ygrad);
   %imshow(phi)
   figure(1)
   %imshow(phi)
   title('phi')
   %figure(2)
   %plot(xgrad)
   %title('xgrad')
   lambda = .24 + n/500000;
   disp(n)
end  



dphi = fatwatfuncs.calcdphi(W, F, theta, phi, S_exp,lambda, zeros(size(S_exp)),Mask);