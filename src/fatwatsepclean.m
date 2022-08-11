
%dbstop if naninf
% load in phantom data and drawn mask
load('spDixonDataFWP_T1_FA10.mat');
load('maskdraw.mat');
rng(1);
S_exp = img(:,:,19);


%initialize parameters
% lambda is this large so that it is on the same scale as phi etc...
lambda =180000;
stepsize = 1;
%theta is given
theta = .2419;

% mask has a threshold and a drawn aspect
Mask = maskdraw .*double(fatwatfuncs.thresholdmasking(abs(S_exp),2000))


%Fat and water are initialized - this could be anything non-zero
% if this is calculated from S_exp, phi is zeroed out.
F = ones(size(S_exp))*15000;
W = ones(size(S_exp))*15000;


% phi is initialized as a guassian filter blur, but could be initialized as
% anything
%phi = linspace(.1,0,192^2)
%phi = reshape(phi,192,192);
phi = imgaussfilt(angle(S_exp),10);
phiinitial = phi;
imshow(phi)
%phi = zeros(size(S_exp))
%dphi = zeros(size(S_exp));




% this just runs until I stop it
for i = 1:500000
    % if phi needs more iterations to converge before W,F is updated
    for n =1:1
        % calculated dphi and update phi
        dphi = fatwatfuncs.calcdphi(W, F, theta, phi, S_exp,lambda, zeros(size(S_exp)),Mask);
        dphi = fatwatfuncs.calcdphi(W, F, theta, phi, S_exp,lambda, dphi,Mask);
        phi = phi + (stepsize)*dphi;
   
        % cost function
        funcs = fatwatfuncs.Psi(W,F,theta,phi,S_exp,lambda,Mask);
      
        figure (1)
        imshow(phi);
    end
%update water and fat
S = S_exp.*exp(-1i*phi);
F = imag(S)/sin(theta);
W = real(S) -F*cos(theta);
end



