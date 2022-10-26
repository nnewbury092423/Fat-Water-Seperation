load('C:\Users\nnewb\Documents\matlab stuff\parameterconfigurations.mat')
filelist = ["ACL11T1_FA15.mat"] %, "ACL12T1_FA15.mat","ACL12T1_FA30.mat"];
phi = [];
x = [];
y = [];
z = [];
masks = cell(1,3);
phases = cell(1,3);
S_exps = cell(1,3) 


for k = 1:numel(filelist)
    %load in the files and configuration
    filename = filelist(k);
    load(filename);
    load('C:\Users\nnewb\Documents\matlab stuff\parameterconfigurations.mat')
    this_record = config(strcmp(filename, {config.key}))
    S_exp = eval(this_record.varname);
    S_exps{k} = eval(this_record.varname);
    [xT,yT,zT] = meshgrid(1:size(S_exps{k},1),1:size(S_exps{k},2),1:size(S_exps{k},3));
    % unwrap phase and create mask
    phase = unwrapPhase(abs(S_exps{k}),angle(S_exps{k}), size(S_exps{k}));
    se = strel('sphere',double(int8(numel(S_exps{k})/80000)));
    Ie = imerode(abs(S_exps{k}),se);
    Iobr = imreconstruct(Ie,abs(S_exps{k}));
    se = strel('sphere',double(int8(numel(S_exps{k})/80000)))
    Iobrd = imdilate(Iobr,se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    Mask1 = (imregionalmax(Iobrcbr)).* double(imfill((abs(S_exp)>7000),"holes"));

%     myfun = 'y~ b0 + b1*x1 + b2*x2 + b3*x3 + b4*x1*x2 + b5*x1*x3 + b6*x2*x3 + b7*(x1^2) +b8*(x2^2)+ b9*(x3^2) + b10*x1*x2*x3 + b11*x1*(x2^2) + b12*(x1^2)*x2 + b13*(x1^3) + b14*(x2^3) + b15*(x3^3)'
%     beta0 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];


    %xyz volume
   

    x = [x; xT(logical(Mask1))];
    y = [y;yT(logical(Mask1))];
    z = [z;zT(logical(Mask1))];

    %phase points
    phi =[phi; phase(logical(Mask1))];
    masks{k} = Mask1
    phases{k} = phase
    % repmat command here
end

myfun = 'y~ b0 + b1*x1 + b2*x2 + b3*x3 + b4*x1*x2 + b5*x1*x3 + b6*x2*x3 + b7*(x1^2) +b8*(x2^2)+ b9*(x3^2) + b10*x1*x2*x3 + b11*x1*(x2^2) + b12*(x1^2)*x2 + b13*(x1^3) + b14*(x2^3) + b15*(x3^3)' %+b16*(x1^4) + b17*(x2^4) + b18*((x1*x2)^2)'
beta0 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];


phiT = reshape(phi,numel(phi),1);
X = [reshape(x,numel(phi),1),reshape(y,numel(phi),1),reshape(z,numel(phi),1)];

mdl = fitnlm(X,phi,myfun,beta0);

XT = [reshape(xT,numel(S_exps{1}),1),reshape(yT,numel(S_exps{1}),1),reshape(zT,numel(S_exps{1}),1)];

phasestuff = reshape(predict(mdl,XT),size(S_exps{1}));
  
M = 0;
S = [];
F = [];
W =[];

for k = 1:numel(filelist)

phasefinal = cell2mat(phases(k)) - phasestuff;
this_record = config(strcmp(filename, {config.key}))
% calculate watind and fatind assumes lowest pixels are water and highest
% pixels are fat
[phaseoff, theta]  = fatwatfuncs.calcoffsettheta((phasefinal.*cell2mat(masks(k))),this_record.thetapercentile,this_record.offsetpercentile,M);


S = cat(4,S,abs(cell2mat(S_exps(k))).*exp(1i*phaseoff));
F = cat(4,F,imag(S(:,:,:,k))/sin(theta));
W = cat(4,W,real(S(:,:,:,k)) -F(:,:,:,k)*cos(theta));

end




