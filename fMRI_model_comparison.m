%% Log evidence model comparison
clear
clc
addpath(genpath('spm12'))
addpath('jbtools')
addpath(genpath('VBA-toolbox-master'))
%% define mask
maskdir ='mask/';
for k = 1:3
    switch k
        case 1
            ROImaskdir = [maskdir,'ACCmask/'];
        case 2
            ROImaskdir = [maskdir,'AICmask/'];
        case 3
            ROImaskdir = [maskdir,'SPLmask/'];
    end
    for s = 1:20
        ROImask(k,s,:) =  scan_nifti_load([ROImaskdir,'subject_0',num2str(s),'.nii']);
    end
end
%% load log evidence
model = {'model(DV)','model(R)','model(DVxCong)','model(DVxconflict1)','model(DVxconflict2)'};
originalpath = 'data/fMRI/';
for k = 1:3  
    k
    for j  = 1:length(model) 
        modelpath =  [originalpath,model{j},'\'];
        for s = 4:23                   
            filename = [modelpath,'subject_0',num2str(s),'_LogEV.nii'];
            templogev = scan_nifti_load(filename);
            LogEvMat(k,s-3,j) = nansum(templogev(logical(ROImask(k,s-3,:))));
        end
    end
end
%% model expected frequencies and exceedance probabilty with VBA toolbox
ROI_name = {'dACC','AIC','SPL'};
for k = 1:3
    [posterior(k),out(k)] = VBA_groupBMC(squeeze(LogEvMat(k,:,:))');
    title(ROI_name{k})
    xticks(1:5)
end