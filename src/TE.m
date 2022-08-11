% Consider multiple peaks for fat
fp = [434 332 486 -76 52];
ap = [0.69 0.13 0.09 0.05 0.04];
TEs = 100
%us
clear('angs','abss')

thetas = 2*pi* (fp) .*(TEs)*1e-6;
disp(thetas)
vecFs = ap.*exp(1j*thetas);
angs =angle(sum(vecFs));
abss =abs(sum(vecFs));
figure(10), subplot(121); 
plot(TEs,angs*180/pi);
subplot(122);
plot(TEs,abss); 
drawnow
figure, plot(TEs,abss); 
ylabel('Signal Intensity (a.u.)'); 
xlabel('TE(us)'); drawnow