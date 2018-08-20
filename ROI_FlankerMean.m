clear 
clc
addpath('function');
addpath('jbtools');
addpath(genpath('spm12'));
%% get contrast
maskdir ='mask/';
path    = 'data/fMRI/GLM4(FM)/';
for k = 1:3
    progressbar(3,k);
    switch k
        case 1
            ROImaskdir = [maskdir,'ACCmask/'];
        case 2
            ROImaskdir = [maskdir,'AICmask/'];
        case 3
            ROImaskdir = [maskdir,'SPLmask/'];
    end
    for s = 1:20
        ROImask =  scan_nifti_load([ROImaskdir,'subject_0',num2str(s),'.nii']);
        tempFM  =  scan_nifti_load([path,'con_subject_0',num2str(s),'.nii']);
        subFM(s,k)         = nanmean(tempFM(logical(ROImask))); 
    end
end
f=figure;
set(f,'Units','inches','position',[0,0,3.5,3]);
barsem(subFM);
xticks(1:3)
xticklabels({'dACC','','AIC','SPL'});