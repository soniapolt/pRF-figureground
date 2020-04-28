function analysisMaster_1_inputFuncs(expt,subjNum)
% input non-retino masks in a flexible way

load(['exptParams/' expt]);
subj = exptSubjs{subjNum};
numFuncs = subjFuncs{subjNum};

funcDir = [fMRIdir expt '/' subj '/funcs'];
outputDir = [fMRIdir expt '/' subj '/matlabAnalysis/data'];
if ~exist(outputDir) mkdir(outputDir); end

%load functional data
for i=numFuncs
    if i <10 ztxt = '0' else ztxt = ''; end
    funcName = [ funcDir '/run' ztxt ' num2str(i) '_filt.nii'];
    
    if ~exist(funcName)
        fprintf('Missing Nifti File %s!\n',funcName);
        break
    end
    [Fdata fhdr] = cbiReadNifti(funcName);
    filename = [outputDir '/funcRun_' num2str(i) '.mat'];
    save(filename,  'fhdr', 'volDimension','Fdata');  
    fprintf('\n\nSaving Func: %s  \n', filename );
end
