
function [vol,siz] = scan_nifti_load(file,mask)
    %% [vol,siz] = SCAN_NIFTI_LOAD(file[,mask])
    % load volumes
    % file : a string, or a cell of strings
    % mask : an array
    % vol  : a volume, or a cell of volumes (vector shaped)
    % siz  : a matrix with the shape of the volumes, or cell of matrices
    % to list main functions, try
    %   >> help scan;
    
    %% function
    
    % mask
    func_default('mask',[]);
    
    % load single file
    if ischar(file)
        vol = spm_vol(file);
        siz = size(vol.private.dat);
        vol = double(vol.private.dat(:));
        if ~isempty(mask)
            assertSize(vol,mask);
            vol(~mask(:)) = [];
        end
        
    % load many volumes
    else
        assert(iscellstr(file),'scan_nifti_load: error. file must be a string or a cell of strings');
        
        % efficient ROI extraction (but doesn't work)
        % if ~isempty(mask)
        %     [x,y,z] = ind2sub(size(mask),find(mask>0));
        %     vol = spm_get_data(file,[x,y,z]')';
        %     vol = mat2vec(num2cell(vol,1));
        %     siz = repmat({[nan,nan,nan]},size(vol));
        %     return
        % end
        
        % recursive call
        vol = cell(size(file));
        siz = cell(size(file));
        for i = 1:numel(file)
            [vol{i},siz{i}] = scan_nifti_load(file{i},mask);
        end
    end

end
