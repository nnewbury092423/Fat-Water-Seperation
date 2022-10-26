   clear
   clc 
   close all


%load in the files and configuration
    filename = 'spDixonDataFWP_T1_FA10.mat';
    load(filename);
    load('parameterconfigurations.mat')
    this_record = config(strcmp(filename, {config.key}))
    S_exp = eval(this_record.varname);
    [xT,yT,zT] = meshgrid(1:size(S_exp ,1),1:size(S_exp,2),1:size(S_exp,3));
    % unwrap phase and create mask
    phase = unwrapPhase(abs(S_exp),angle(S_exp), size(S_exp));
    Mask = roipoly(abs(S_exp(:,:,25)));
    Mask1 = zeros(size(S_exp));
    Mask1(:,:,19:25) = repmat(Mask,1,1,7);
    Mask1(mean(abs(S_exp) + 5> abs(S_exp) > mean(abs(S_exp) - 5)))

    myfun = 'y~ b0 + b1*x1 + b2*x2 + b3*x3'
    beta0 = [0,0,0,0];


    %xyz volume
   

    x = xT(Mask1>0);
    y = yT(Mask1>0);
    z = zT(Mask1>0);

    %phase points
    phi=phase(Mask1>0);
    %masks{k} = Mask
    %phases{k} = phase
    % repmat command here

phiT = reshape(phi,numel(phi),1);
X = [reshape(x,numel(phi),1),reshape(y,numel(phi),1),reshape(z,numel(phi),1)];

mdl = fitnlm(X,phi,myfun,beta0);

XT = [reshape(xT,numel(S_exp),1),reshape(yT,numel(S_exp),1),reshape(zT,numel(S_exp),1)];

phasestuff = reshape(predict(mdl,XT),size(S_exp)).*Mask1;


phasefinal = (phase - phasestuff).*Mask1;

noiseelem = phasefinal(abs(phasefinal)>0)

histogram(noiseelem,500);

sigma = std(noiseelem);
mean = mu(noiselem);
mag = mean(abs(S_exp(Mask1>0)));
Title(['Histogram of Phase for Homogenous Area Signal Mag = ' mag])
xlabel('Phase Mu = 0')
ylabel('Frequency')