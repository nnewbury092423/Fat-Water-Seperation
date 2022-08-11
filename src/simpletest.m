x = 4*exp(1i*-.1);
A = 4*exp(1i*.3);
phi = 0;
%A = abs(Sexp).*exp(1i*theta)
%initialize
%F = [1 3; 2 4]; %abs(Sexp) %imag(Sexp)/sin(theta);
%W = [3 4.1; 3.2 2.5];   %real(Sexp) - abs(F)*cos(theta);
%F = ones(3,3)
%W = ones(3,3)
%A = W + F*exp(1i*theta);
%amp = abs(A);
%aphase = angle(A);
%theta = .326;
%actphase = [.12, .35, .23 ;.3, .48, .33;.67, 1.5, .76];
%S_exp = [ 4 5 8.9;1 2 3.1;12 6 5.4].*exp(1i*actphase);
phi = 0;
for n =1:100
    dphi = -imag(A*exp(1i*phi)*conj(x))/(abs(A)^2);
    phi = phi + dphi;
    disp(phi)
    %S = S_exp.*exp(-1i*phi)
    %F = imag(S)/sin(theta);
    %W = real(S) - abs(F)*cos(theta);
    watfatphase = angle(x) - phi;
    fatwatfuncs.Psiamp(A, phi,x);
end  