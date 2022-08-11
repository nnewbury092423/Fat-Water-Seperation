clear
close all
clc
rng(1)
sizenum = 100
Wat=abs(randn(sizenum)*1000);
Fat=abs(randn(sizenum)*1000);
theta = .2419
% phi is a gradually changing error
phisig = linspace(.5,-.5,sizenum^2);
phisig = reshape(phisig,sizenum,sizenum);
stepsize = 1;

S_exp = (Wat + Fat*exp(1i*theta)).*exp(1i*phisig);
dcouplesignal = (Wat + Fat*exp(1i*theta));

phi = imgaussfilt(angle(S_exp),10);
lambda = 110
%S = S_exp.*exp(-1i*phi)*exp(1i*.121);
F = ones(size(S_exp))*30;%imag(S)/sin(theta);
W = ones(size(S_exp))*30;%real(S) - F*cos(theta);
F = abs(imag(S_exp)/sin(theta));
W = abs(real(S_exp) - F*cos(theta));
Mask =ones(size(S_exp));%double(fatwatfuncs.thresholdmasking(S_exp,5000));
disp(phi(9,76))
  for i = 1:500000
    for n =1:10
        disp('phi pre')
        disp(phi(9,76))
        disp('what it is supposed to be similar to')
        disp(angle(S_exp(9,76)))
        S_temp = (W + F*exp(1i*theta)).*exp(1i*phi);
         dphi = fatwatfuncs.calcdphi(W, F, theta, phi, S_exp,lambda, zeros(size(S_exp)),Mask);
        [dphi, phigradient, dphigradient]  = fatwatfuncs.calcdphi(W, F, theta, phi, S_exp,lambda, dphi,Mask);
        phi = phi +(stepsize)*dphi;
   
        funcs = fatwatfuncs.Psi(W,F,theta,phi,S_exp,lambda,Mask);
       
        %S = S_exp.*exp(-1i*phi);
        %F = imag(S)/sin(theta);
        %W = real(S) - F*cos(theta);
        figure (1)
        imshow(phi);
        disp('phi post')
        disp(phigradient(9,76))
        disp(dphigradient(9,76))
        disp(phi(9,76))
        
        %disp('here')
    end
S = S_exp.*exp(-1i*phi);
F = abs(imag(S)/sin(theta));
W = abs(real(S) -F*cos(theta));
%disp('down')
disp(phi(80,80))
%disp('same')
%disp(phi(147,60))
end


s=diag(S);
A=U*diag(s+max(s))*U'; % to make A symmetric, well-contioned
b=randn(1000,1);

x=conjgrad(S_exp,vals)


tic,x1=A\b;toc

norm(x-x1)

norm(x-A*b)