clear
clc
addpath(genpath('spm12'));
addpath('jbtools');
addpath('function');
%% define mask
maskdir ='mask/';
filedir ='data/fMRI/GLM(binned)/';
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
    for t = 1:4
        for f = 1:4
            foldname =['Tm',num2str(t),'_Fm',num2str(f),'_Fv/'];
            fulldir  =[filedir,foldname];
            for s = 4:23                
                if s <10; substr = 'subject_00';
                else                    
                    substr = 'subject_0';
                end
                ROImask  =  scan_nifti_load([ROImaskdir,'subject_0',num2str(s-3),'.nii']);
                Files = dir([fulldir,substr,num2str(s),'*.nii']);
                allvx    = scan_nifti_load([fulldir,Files.name]);
                subTM(k,s-3,t,f)  = nanmean(allvx(logical(ROImask)));
            end
        end
    end
end
%% plot
ROI_name = {'dACC','AIC','SPL'};
f = figure;
set(f,'Units','inches','position',[0,0,10,3])
for k = 1:3
    subplot(1,3,k);
    linesem(squeeze(subTM(k,:,:,:)));
    xlim([0.5 4.5])
    title(ROI_name{k})
    xlim([0.5 4.5]);
    xticks([1:4])
    xticklabels({'-High','-Low','+Low','+High'});
    xlabel('Flanker mean orientation')    
    ylabel('parameter estimates')
    
end