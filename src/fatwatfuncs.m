classdef fatwatfuncs
    properties
        Mask
    end
    methods(Static)
        % returns Laplacian of 2 dimensional image. Used to take Laplacian
        % of phi and dphi
        function result = Smoothing2d(phi)

            xgrad = [zeros(length(phi(1,:)),1) diff(phi,1,2)];
            ygrad = [zeros(1, length(phi(:,1))); diff(phi,1,1)];
            %disp(mean((xgrad+ygrad),'all'))
            xxgrad = [diff(xgrad,1,2) zeros(length(phi(1,:)),1)];
            yygrad = [diff(ygrad,1,1); zeros(1,length(phi(:,1)))];
            %[yxgrad, yygrad] = gradient(ygrad);
            result = (xxgrad + yygrad);
        end

        function result = Smoothing3d(phi)
            %xgrad = [zeros(length(phi(1,:,:)),1) diff(phi,1,2)];
            
            diffx = diff(phi,1,2);
            zerosx = zeros(size(phi(:,1,:)));
            xgrad = cat(2,zerosx,diffx);

            diffy = diff(phi,1,1);
            zerosy = zeros(size(phi(1,:,:)));
            ygrad = cat(1,zerosy,diffy);

            diffz = diff(phi,1,3);
            zerosz = zeros(size(phi(:,:,1)));
            zgrad = cat(3,zerosz,diffz);
            %disp(mean((xgrad+ygrad),'all'))
            
            diffxx = diff(xgrad,1,2);
            zerosxx = zeros(size(xgrad(:,1,:)));
            xxgrad = cat(2,diffxx,zerosxx);

            diffyy = diff(ygrad,1,1);
            zerosyy = zeros(size(ygrad(1,:,:)));
            yygrad = cat(1,diffyy,zerosyy);

            diffzz = diff(zgrad,1,3);
            zeroszz = zeros(size(zgrad(:,:,1)));
            zzgrad = cat(3,diffzz,zeroszz);




            result = (xxgrad + yygrad + zzgrad);
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
                sum(max(0,-(angle(S_exp)-phi)).^2,'all') + sum(max(0,angle(S_exp)-phi-theta).^2,'all');
        end

        % calculates dphi g
        function [dphi, phigradient,dphigradient] = calcdphi(W,F,theta,phi,S_exp,lambda,dphi,Mask,dim)
            if  dim == 2
                phigradient = 2*lambda*fatwatfuncs.Smoothing2d(phi);
                dphigradient = 2*lambda*fatwatfuncs.Smoothing2d(dphi);
            elseif dim == 3
                 phigradient = 2*lambda*fatwatfuncs.Smoothing3d(phi);
                dphigradient = 2*lambda*fatwatfuncs.Smoothing3d(dphi);
            end

            mag = (W.^2 + F.^2 + 2*W.*F*cos(theta));
            phi1gradient = ones(size(S_exp)).* double(angle(S_exp)-phi-theta>0)*(4*lambda).*-(angle(S_exp)-phi-theta);
            phi2gradient = ones(size(S_exp)).* double(angle(S_exp)-phi<0)*(4*lambda).*(angle(S_exp)-phi);


            imagpart = imag((W + F*exp(1i*theta)).*exp(1i*phi).*conj(S_exp));
            %disp(mean(imagpart, 'all'))

            dphi = (-imagpart + Mask.*(phigradient + dphigradient + phi1gradient+ phi2gradient))./(mag);
        end

        function Mask = thresholdmasking(S,T)
            Mask = (S < T);
        end
    end
end