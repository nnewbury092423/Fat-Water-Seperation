load('spDixonDataFWP_T1_FA10.mat');
load('maskdraw.mat');
hold off
Mask = double(maskdraw)
S_exp = double(img(:,:,19));
arrlength = length(S_exp(:,1));
lambda = 1;
theta = .2419;
n = 200;
cost = zeros(n,1); 
% make a phi array
%flesh out water and fat

phioff = linspace(0,1,n);

%3 dimensional matrix
phi = ones(size(S_exp));

phiarr = ones(arrlength ,arrlength , n);

count = 1;
minw  = zeros(1,n)
minf = zeros(1,n)
for phi = linspace(-1,1,n)
       
    phiarr(:,:,count) = ones(size(S_exp)) * phi;
    
    S = S_exp*exp(-1i*phi);
    F = imag(S)/sin(theta);
    W = real(S) - F*cos(theta);
    minf(count) = min(min(F));
    minw(count) = min(min(W));
    cost(count) = fatwatfuncs.Psitot(W,F,theta,phiarr(:,:,count),S_exp,lambda,Mask);
  
    count  = count+1;
end
phi = reshape(phiarr(1,3,:),1,n)
 figure(1) 
 plot(phi,cost)

 
 figure(2)
 hold on
 plot(phi, minw)
 plot(phi,minf)
 
 % limit water and fat to negative to positive
 
 % graph showing phi and total cost

 % phi map and cost map