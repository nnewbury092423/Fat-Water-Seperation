classdef fatwatfuncs
    properties
        Mask
    end
    methods(Static)
          function result = Smoothing(phi)
          xgrad = [zeros(length(phi(1,:)),1) diff(phi,1,2)];
          ygrad = [zeros(1, length(phi(:,1))); diff(phi,1,1)];
          %disp(mean((xgrad+ygrad),'all'))
          xxgrad = [diff(xgrad,1,2) zeros(length(phi(1,:)),1)];
          yygrad = [diff(ygrad,1,1); zeros(1,length(phi(:,1)))];
          %[yxgrad, yygrad] = gradient(ygrad);   
          result = (xxgrad + yygrad);
          
           
          end
          function result = Psi(W,F,theta,phi,S, lambda,Mask)
              [xgrad, ygrad] = gradient(phi);
              result = .5 * abs((W + F.*exp(1i*theta)).*exp(1i*phi) - S).^2 ...
              + lambda*Mask.*( xgrad^2 + ygrad^2);          
          end
          function result = Psiamp(A, phi,Sexp)
              result= .5 *abs((A.*exp(1i*phi) - Sexp)).^2;
              
          end
          function dphi = calcdphi(W,F,theta,phi,S_exp,lambda,dphi,Mask)
              mag = (W.^2 + F.^2 + 2*W.*F*cos(theta));
              
              phigradient = 2*lambda*fatwatfuncs.Smoothing(phi);
              dphigradient = 2*lambda*fatwatfuncs.Smoothing(dphi);
              %disp('phigradient')
             % disp(mean(phigradient,'all'))
              imagpart = imag((W + F*exp(1i*theta)).*exp(1i*phi).*conj(S_exp));
              %disp(mean(imagpart, 'all'))
             
              dphi = (-imagpart + Mask.*(phigradient + dphigradient))./(mag);% ...
          end
          
          function Mask = thresholdmasking(S,T)
                Mask = (S > T);
          end
    end
end