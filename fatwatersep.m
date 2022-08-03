
%dbstop if naninf
% load data in
load('spDixonDataFWP_T1_FA10.mat');
load('maskdraw.mat');
rng(1);
S_exp = img(:,:,19);

% guassian filter
%initialize parameters
lambda =150000;
stepsize = 1;
theta = .2419;
%Mask = double(fatwatfuncs.thresholdmasking(S_exp,5000));
Mask = maskdraw .*double(fatwatfuncs.thresholdmasking(S_exp,1000))
%initialize phi, dphi with fat and water guesses
%F =  double(imag(S_e n xp)/sin(theta)>0).*imag(S_exp)/sin(theta) +1;
%W =  double((real(S_exp)- F*cos(theta))>0).*(real(S_exp) - F*cos(theta)) + 1;

F = ones(size(S_exp))*15000;%imag(S_exp)/sin(theta);
W = ones(size(S_exp))*15000;%real(S_exp) - F*cos(theta);

%phi = rand(192);
%dphi =rand(192);
%phi = linspace(.1,0,192^2)
%phi = reshape(phi,192,192);
phi = imgaussfilt(angle(S_exp),10);
phiinitial = phi;
imshow(phi)
% %dphi = zeros(size(S_exp));
%val = (W + F*exp(1i*theta));
%ang = angle(val);
%S = S_exp.*exp(-1i*phi);
%F = imag(S)/sin(theta);
%W = real(S) - F*cos(theta);
for i = 1:500000
    for n =1:1
        dphi = fatwatfuncs.calcdphi(W, F, theta, phi, S_exp,lambda, zeros(size(S_exp)),Mask);
        dphi = fatwatfuncs.calcdphi(W, F, theta, phi, S_exp,lambda, dphi,Mask);
        phi = phi + (stepsize)*dphi;

        funcs = fatwatfuncs.Psi(W,F,theta,phi,S_exp,lambda,Mask);
       
        %S = S_exp.*exp(-1i*phi);
        %F = imag(S)/sin(theta);
        %W = real(S) - F*cos(theta);
        figure (1)
        imshow(phi);
        %disp('here')
    end
S = S_exp.*exp(-1i*phi);
F = imag(S)/sin(theta);
W = real(S) - F*cos(theta);
disp('down')
disp(phi(120,80))
disp('same')
disp(phi(147,60))
end

%S = S_exp.*exp(-1i*phi);
%F = imag(S)/sin(theta);
%W = real(S) - abs(F)*cos(theta);


% initialize fat and water using initial phi
% S = S_exp.*exp(-1i*phi);
% F_exp = imag(S)/sin(theta);
% W = real(S) - abs(F)*cos(theta);
% %initialize phi, dphi with fat and water guesses
%  for n =1:200
%      dphi = fatwatfuncs.calcdphi(W, F, theta, phi, S_exp, lambda,dphi,Mask);
%      phi = phi + (stepsize)*dphi;
%      imshow(phi)
%      S = S_exp.*exp(-1i*phi);
%      F = imag(S)/sin(theta);
%      W = real(S) - abs(F)*cos(theta);
%  end
 %phi1 = filter(100, 1, phi)


% %phi = zeros(size(S_exp));
% %dphi = zeros(size(S_exp));
% for n = 1:100
%    %update dphi
%    phi = zeros(size(S_exp));
%    dphi = zeros(size(S_exp));
%     for n =1:20
%         dphi = fatwatfuncs.calcdphi(W, F, theta, phi, S_exp, lambda,dphi,Mask);
%         phi = phi + (stepsize)*dphi;
%     end
%    %update F,W
%    
%     S = S_exp.*exp(-1i*phi);
%     F = imag(S)/sin(theta);
%     W = real(S) - abs(F)*cos(theta);
%     funcs = fatwatfuncs.Psi(W,F,theta,phi,S_exp,lambda,Mask);
%     %imshow(phi) 
% end
% 
% phi2 = (phi1 +phi)/2;
% 
% 


