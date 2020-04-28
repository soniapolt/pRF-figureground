
clear all; close all;

expt = 'LGNedgeAll'; % experimental details for main expt + gap control (collected in same session)

load(['exptParams/' expt]);

% which subjects are we including?

subjNums = [1:length(exptSubjs)]; % all of them!

labelNames = {'lh_LGN' 'rh_LGN'}; % loaded in as mask.nii files, since they are defined in the volume not on surface

hems = {'lh' 'rh'};
ROIs =  {'LGN' 'V1' 'V2' 'V3'}

percCutoff = [67 67 67 67]; % r2 percentile cutoff for pRFs in each ROI

for s = subjNums
    
    %%% data loading: relies on cbiReadNifti.m from http://gru.stanford.edu/mrTools
    % input LGN masks
    analysisMaster_1_inputMasks(expt,s,labelNames);
    % input functional data
    analysisMaster_1_inputFuncs(expt,s);
    
    for h = 1:length(hems)
        hem = hems{h};
        for r = 1:length(ROIs)
            ROI = ROIs{r};
            if strcmp(ROI,'LGN')==0
                maskPre = 'all_'; else maskPre = ''; end
            
            % read the voxelwsie data (no ROI averaging)
            analysisMaster_2_voxelWiseRead(expt,s,hem,ROI,maskPre)
            
            % run a GLM analysis across all runs
            voxelwiseMaster_3_GLM_blockConds(expt,s,hem,ROI)
        end
    end
    % make ROIs by whether pRF is in/out of stim
    % grabs percentile of each subject's voxels in
    % each ROI
    prfMaster_4_GLM_sortVox(expt,s,ROIs,percCutoff)
end

prfMaster_5_GLM_groupSortVox(expt,ROIs,percCutoff)

cutoffs = 67;
% %%% use ridge regression to estimate the spatial profile of effects
for r= 1:length(ROIs)
    for c = cutoffs
        voxelwiseMaster_4_subjRidge(expt,subjNums,ROIs{r},c,['cvRidge_' num2str(c)],0);
    end
end

for r = 1:length(ROIs)
    prfMaster_5_GLM_groupInOut_ROIacrossSubjs(expt,ROIs{r},percCutoff(r))
end
