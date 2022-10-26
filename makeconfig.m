
fnames = {'key','varname', 'se_ratio', 'thetapercentile','offsetpercentile','watind','fatind', 'mag'};
content = {'T1_FA15.mat','img_Full_c',80000, 10,80, [99,70,25], [120,82,25],4000;  ...
           'T1_FA5.mat','img_Full_c',80000, 10,80, [99,70,25], [120,82,25],4000; ...
            'P10752.7_recon_nom_1pD.mat_SP-Dixon_knee.mat','S11',80000, 7,50, [99,70,25], [120,82,25],4000; ...
            'P44544.7_recon_nom_1pD.mat_SP-Dixon_ankle.mat','S11',80000, 10,80, [99,70,25], [120,82,25],4000; ...
            'spDixonDataFWP_T1_FA10.mat','img',80000, 10,80,[99,70,25], [120,82,25],4000; ...
            'ACL11T1_FA5.mat','img_Full_c', 80000, 20,50, [99,70,25], [120,82,25],4000;  ...
            'ACL11T1_FA15.mat','img_Full_c',80000, 20,50, [99,70,25], [120,82,25],4000;  ...
            'ACL11T1_FA30.mat','img_Full_c',80000, 20,50, [99,70,25], [120,82,25],4000; ...
                    'ACL12T1_FA5.mat','img_Full_c', 80000, 20,50, [99,70,25], [120,82,25],4000;  ...
            'ACL12T1_FA15.mat','img_Full_c',80000, 20,50, [99,70,25], [120,82,25],4000;  ...
            'ACL12T1_FA30.mat','img_Full_c',80000, 20,50, [99,70,25], [120,82,25],4000};
config = cell2struct(content, fnames, 2);

save('parameterconfigurations', 'config')
disp('here')