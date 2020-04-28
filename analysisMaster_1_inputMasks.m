function analysisMaster_1_inputMasks(expt,subjNum,labelNames)
% input non-retino masks in a flexible way

load(['exptParams/' expt]);

subj = exptSubjs{subjNum};

maskDir = [fMRIdir  expt '/' subj '/masks'];

%load masks/labels

for i=1:length(labelNames)
    maskName = [ maskDir '/' labelNames{i} '.mask.nii'];
    [Fdata fhdr] = cbiReadNifti(maskName);
    tmpdata = Fdata(:)>0; % 0,1,2,3 --> ones and zeros
    vox_ind = find(tmpdata==1);
    nVox = length(vox_ind);
    
    Format = fhdr.matlab_datatype; % int32
    
    data = zeros(size(Fdata)); % 64   28    64
    data(vox_ind) = 1;
    mask = data;
    maskData = [maskDir '/' labelNames{i} '.mat'];
    fprintf('\n\nSaving Mask: %s  \n', maskData );
    save(maskData, 'mask');
    
%     scaling data from double to int32
%    data = eval( [ Format '(data)' ] ); 
end
