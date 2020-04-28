function analysisMaster_2_voxelWiseRead(expt,subjNum,hem,ROI,maskPre)
% this analysis calculated the mean PSC relative to each run, not the whole
% experiment. it also allows us to plot time series.
% 8/18/14 update: for objectPRF fitting
% PSC,  which we can then shift relative to the baseline of blanks
% and raw data renormalized to baseline of blanks (like old prewindow norm)
%clear all; %close all;

load([pwd '/exptParams/' expt '.mat']);

subj = exptSubjs{subjNum};
numFuncs = subjFuncs{subjNum};

analysis.whichMask = maskPre;
analysis.trimOutsideVol =1; % trim voxels that leave the volume
analysis.trimOutliers = 1; % trim voxels that fluctuate above a certain threshold
analysis.outlierCutoff = 3; %in SD's, for each voxel

% visualize runwise mean BOLD as a basic QA in each subject
if strcmp(ROI,'V1') plotRunMeans = 1; else plotRunMeans = 0; end

% before you start, copy the functional runs' matlab output from the scan
% into a new 'scannerOutput' directory, and rename them as run##_scannerOutput.mat

dataDir = [fMRIdir expt '/' subj '/matlabAnalysis/data'];
scanOutputDir = [fMRIdir expt '/' subj '/scan_matlabOutput'];
maskDir = [fMRIdir expt '/' subj '/masks'];
outputDir = [fMRIdir expt '/' subj '/matlabAnalysis/voxelwiseData'];
if ~exist(outputDir) mkdir(outputDir); end

% load basic params, which should be the same for both the polar and the eccen runs
% load([scanOutputDir '/run1_scanOutput.mat']);

% functional name convention: funcRun_##.mat
% it contains: Fdata

% mask name convention: {labelNames}.mat
% it contains: mask

%%% load mask
load([maskDir '/' analysis.whichMask  hem ROI '.mat']);

% reshape the mask so we can apply it later
analysis.voxelsInROI = sum(mask(:));
%numVoxels=analysis.volDimension(1)*analysis.volDimension(2)*analysis.volDimension(3);

for a = numFuncs
    % load the data 
    load([dataDir '/funcRun_' num2str(a) '.mat']);
    analysis.exptTRs = size(Fdata,4);
    analysis.volDimension = [size(Fdata,1) size(Fdata,2) size(Fdata,3)];
    
    if analysis.exptTRs < 20
       fprintf('Uh oh! Data is too small to make sense! Quitting out now. \n');
       break
    end
    
    % if this code does locs & runs, which are of variable TR lengths, we
    % need to do this step in the for loop now...
    analysis.voxelInd = []; reshapedMask = []; maskPattern = [];
    maskPattern  = repmat(mask,[1,1,1,analysis.exptTRs]); % this is volDims(1),volDims(2),volDims(3),numTRs size, just like our volume

    reshapedMask = reshape(maskPattern,prod(analysis.volDimension),analysis.exptTRs);
    [analysis.voxelInd,~] = find(reshapedMask(:,1)>0);

    fprintf('Working on Run %d...\n',a)
    voxelWise=struct('volumeInd',[],'runNum',a,'PSC',{zeros(1,analysis.exptTRs)});
    
    for n = 1:analysis.voxelsInROI
        voxelWise(n).volumeInd = analysis.voxelInd(n);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % apply mask and reshape data                                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % collapse the data into voxels x TR
    reshapedData = reshape(Fdata,prod(analysis.volDimension),analysis.exptTRs);
    [analysis.voxelInd,time] = find(reshapedMask(:,1)>0);
    reshapedData= reshapedData(analysis.voxelInd,:);% trim voxels that have been zero'ed out by our mask
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % trim voxels that leave the volume                                                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if analysis.trimOutsideVol ==1
        % cut out zero voxels
        [analysis.zeroVoxels,time] = find(reshapedData<=0);
        analysis.zeroVoxels = unique(analysis.zeroVoxels);
        for m = 1:length(analysis.zeroVoxels)
            reshapedData(analysis.zeroVoxels(m),:)= nan(1,analysis.exptTRs);
        end
        
    end
    analysis.percentOutsideVol = 100*length(analysis.zeroVoxels)/analysis.voxelsInROI;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % trim outlier timepoints at each voxel                                                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if analysis.trimOutliers == 1
        % calculate mean and std of every voxel in the ROI
        meanOfVoxels = nanmean(reshapedData,2);
        stdOfVoxels = std(reshapedData,0,2);
        topEachVoxel = meanOfVoxels+analysis.outlierCutoff*stdOfVoxels;
        bottomEachVoxel = meanOfVoxels-analysis.outlierCutoff*stdOfVoxels;
        
        %check for threse 3STD outliers, and replace them with the cutoff
        analysis.allOutliers = [];
        for n=1:length(topEachVoxel)
            for m = 1:analysis.exptTRs
                if reshapedData(n,m)>topEachVoxel(n) 
                    reshapedData(n,m) = topEachVoxel(n);
                elseif reshapedData(n,m)<bottomEachVoxel(n)
                    reshapedData(n,m) = bottomEachVoxel(n);
                    analysis.allOutliers = [analysis.allOutliers;n,m];
                end
            end
        end
        analysis.numOutliers = length(analysis.allOutliers);
        analysis.percentOutlier = 100*analysis.numOutliers/(analysis.voxelsInROI*analysis.exptTRs);
    end
    
    % to calculate the percent signal change (PSC)
    meanData = nanmean(reshapedData,2); % one value for each voxel across TRs
    meanPattern  = repmat(meanData,1,analysis.exptTRs);
    PSC = (reshapedData - meanPattern)./meanPattern; % calculate average percent signal change across the volume for this run
    
    Fdata = [];maskedData = [];meanPattern = [];meanData = [];
    
    % PSC is size numVoxels x numTRs
    % meanPSC is size 1 x numTRs;
    
    for v = 1:analysis.voxelsInROI
      voxelWise(v).PSC = PSC(v,:);
      voxelWise(v).rawData = reshapedData(v,:); 
    end

    fprintf('Saving...\n');
    
    if plotRunMeans ==1
    % plot run mean
    runMean = [];
    for v = 1:analysis.voxelsInROI
       runMean = [runMean; voxelWise(v).PSC]; 
    end
    figure;  plot(nanmean(runMean)); hold on;
       vline([8:16:length(voxelWise(v).PSC)],'r');
       title(['Run ' num2str(a) ' ' ROI]);
    end
    
    % save the voxelWise data
    eval(['save ([outputDir ''/voxelWise_run' num2str(a) '_' hem ROI '.mat'' ],''analysis'',''voxelWise'',''subj'')']); 
end
end

