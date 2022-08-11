classdef fatwatfuncs
    properties
        Mask
    end
    methods(Static)
        % returns Laplacian of 2 dimensional image. Used to take Laplacian
        % of phi and dphi
          function result = Smoothing(phi)
          xgrad = [zeros(length(phi(1,:)),1) diff(phi,1,2)];
          ygrad = [zeros(1, length(phi(:,1))); diff(phi,1,1)];
          %disp(mean((xgrad+ygrad),'all'))
          xxgrad = [diff(xgrad,1,2) zeros(length(phi(1,:)),1)];
          yygrad = [diff(ygrad,1,1); zeros(1,length(phi(:,1)))];
          %[yxgrad, yygrad] = gradient(ygrad);   
          result = (xxgrad + yygrad);
          
           
          end
           function result = Differ(phi)
          xgrad = [zeros(length(phi(1,:)),1) diff(phi,1,2)];
          ygrad = [zeros(1, length(phi(:,1))); diff(phi,1,1)];
          result = (xgrad + ygrad);
          
           
          end
          function result = Psi(W,F,theta,phi,S, lambda,Mask)
              [xgrad, ygrad] = gradient(phi);
              result = .5 * abs((W + F.*exp(1i*theta)).*exp(1i*phi) - S).^2 ...
              + lambda*Mask.*( xgrad^2 + ygrad^2);          
          end
 
          function phi = wrap(phi)
              phi = rem(phi,2*pi);
              
              
              phi = rem(phi - pi,-pi);
              
          end
          function result = Psiamp(A, phi,S_exp)
              [xgrad, ygrad] = gradient(phi);
              result = .5 *abs((A.*exp(1i*phi) - S_exp)).^2 + lambda*Mask.*( xgrad^2 + ygrad^2);
              
          end
          
          function result = Psiamptot(A,phi,S_exp,lambda,Mask)
              [xgrad, ygrad] = gradient(phi);
              result = .5 *sum(abs((A.*exp(1i*phi) - S_exp)).^2,'all') + sum(lambda*Mask.*( xgrad.^2 + ygrad.^2),'all');
              
          end
          
         
          function result = Psitot(W,F,theta,phi,S_exp, lambda,Mask)
              [xgrad, ygrad] = gradient(phi);
              result = .5 * sum(abs((W + F.*exp(1i*theta)).*exp(1i*phi) - S_exp).^2,'all') ...
              + sum(lambda*Mask.*( xgrad^2 + ygrad^2),'all')+ ...
              sum(max(0,-(angle(S_exp)-phi)).^2,'all') + sum(max(0,angle(S_exp)-phi+theta).^2,'all');          
          end
         
          % calculates dphi g
          function [dphi, phigradient,dphigradient] = calcdphi(W,F,theta,phi,S_exp,lambda,dphi,Mask)
              mag = (W.^2 + F.^2 + 2*W.*F*cos(theta));
              
              phi1gradient = ones(size(S_exp)).* double(angle(S_exp)-phi-theta>0)*(1*lambda).*(angle(S_exp)-phi-theta);
              phi2gradient = ones(size(S_exp)).* double(angle(S_exp)-phi<0)*(1*lambda).*(angle(S_exp)-phi);
              
              %               if (-(angle(S_exp)-phi)<0)
%                   phi2gradient = 0;
%               else
%                   phi2gradient = 2;
%               end
%               
              
             % Fgradient
              phigradient = 2*lambda*fatwatfuncs.Smoothing(phi);
              dphigradient = 2*lambda*fatwatfuncs.Smoothing(dphi);
              %disp('phigradient')
             % disp(mean(phigradient,'all'))
              imagpart = imag((W + F*exp(1i*theta)).*exp(1i*phi).*conj(S_exp));
              %disp(mean(imagpart, 'all'))
              
              dphi = (-imagpart + Mask.*(phigradient + dphigradient + phi1gradient+ phi2gradient))./(mag); 
          end
          
          function Mask = thresholdmasking(S,T)
                Mask = (S < T);
          end
    end
end