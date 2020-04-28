function voxelwiseMaster_3_GLM_blockConds(expt,subjNum,hem,ROI,normOff,whichHRF)
% 10/18/16 edit: makes this all very flexible. can, in theory, take any
% experiment with parameters written in exptParams and chop up data by
% conditions, fit block-wise GLM


load([pwd '/exptParams/' expt]);
subj = exptSubjs{subjNum};
numRuns = subjFuncs{subjNum};

scanOutputDir = [fMRIdir '/' expt '/' subj '/scan_matlabOutput'];
dataDir = [fMRIdir '/' expt '/' subj '/matlabAnalysis/voxelwiseData'];
GLMrunDir = [fMRIdir '/' expt '/' subj '/matlabAnalysis/GLMrunData'];
if ~exist(GLMrunDir) mkdir(GLMrunDir); end
GLMresultsDir = [fMRIdir '/' expt '/' subj '/matlabAnalysis/GLMresults'];
if ~exist(GLMresultsDir) mkdir(GLMresultsDir); end

% GLM result struct - across whole experiment
load([dataDir '/voxelWise_run' num2str(numRuns(1)) '_' hem ROI '.mat']);

c = struct('condName',condNames,'blockBetas',[]);
GLM(1:analysis.voxelsInROI) = struct('volumeInd',[],'betas',[],'blanks',[],'cond',c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make timing file for each run, fit GLMs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for currentRun = numRuns
    
    % load the data and timing stuff for this run
    load([dataDir '/voxelWise_run' num2str(currentRun) '_' hem ROI '.mat']);
    
    % timing for this run
    load([scanOutputDir '/runOutput_' num2str(currentRun) '.mat']);
    
    % data struct for this run
    
    % preallocating.
    GLMdata(1:analysis.voxelsInROI) = struct('PSCdata',[],'zData',[],'volumeInd',[],'betas',[],'betaInts',[],'resids',[],'stats',[],'stErrBetas',[]);
    
    % housekeeping - scanner output is in seconds, data is in TRs
    blockStarts = experiment.startBlock/2+1;
    numTRs = length(voxelWise(1).PSC);
    
    blockNum = 0; condOrder = [];
    GLMmodel = zeros(numTRs,length(blockStarts));
    
    for n = 1:length(conditions) % extract congruent blocks
        for m = 1:length(conditions(n).startTimes)
            startTR = conditions(n).startTimes(m)/2+1;
            endTR = startTR + params.blockLength/params.TRlength-1;
            blockNum = blockNum+1;
            
            thisCond = analysisMaster_condAssignment(expt,conditions,experiment,hem,n,m);
            
            GLMmodel(startTR:endTR,blockNum)=1; % fill model - seperate regressor for each block
            condOrder = [condOrder thisCond];
        end
    end
    
    %figure;plot(GLMmodel);
    
    for n=1:blockNum
        switch whichHRF
            case 'spm'
                conved = conv(spm_hrf(2),GLMmodel(:,n));
            case 'my'
                conved = conv(HIRF_doubleGamma,GLMmodel(:,n));
        end
        GLMmodel(:,n)=conved(1:numTRs,:); %convolves each run with the double gamma HRF
    end
    GLMmodel = [GLMmodel ones(numTRs,1)];
    
    for v = 1:analysis.voxelsInROI
        GLMdata(v).volumeInd = voxelWise(v).volumeInd;
        GLMdata(v).PSCdata = [GLMdata(v).PSCdata voxelWise(v).PSC];
        
        % zscore calc
        zData  = (voxelWise(v).rawData - mean(voxelWise(v).rawData))./std(voxelWise(v).rawData);
        GLMdata(v).zData = [GLMdata(v).zData zData];
        
        %%% REGRESS
        GLM(v).volumeInd = voxelWise(v).volumeInd;
        
        if ~exist('normOff','var')||normOff == 0
            [GLMdata(v).betas,GLMdata(v).betaInts,GLMdata(v).resids,~,GLMdata(v).stats]=regress(GLMdata(v).zData',GLMmodel);
        else
            [GLMdata(v).betas,GLMdata(v).betaInts,GLMdata(v).resids,~,GLMdata(v).stats]=regress(GLMdata(v).PSCdata',GLMmodel);
        end
        
        % aggregate betas by condition - across runs
        for c = 1:length(condOrder)
            GLM(v).cond(condOrder(c)).blockBetas = [GLM(v).cond(condOrder(c)).blockBetas GLMdata(v).betas(c)];
        end
        
        GLM(v).blanks = [GLM(v).blanks GLMdata(v).betas(end)];
        
    end
    if ~exist('normOff','var')||normOff == 0
        eval(['save ([GLMrunDir ''/voxelWise_GLMdata_run' num2str(currentRun) '_' hem ROI '.mat'' ],''GLMdata'')']);
    else
        eval(['save ([GLMrunDir ''/voxelWise_GLMdata_run' num2str(currentRun) '_' hem ROI '_nsBetas.mat'' ],''GLMdata'')']);
    end
end % all runs

for v = 1:analysis.voxelsInROI
    for c = 1:length(condNames)
        GLM(v).betas(c) = mean(GLM(v).cond(c).blockBetas);
        voxBetas(v,c) = GLM(v).betas(c);
    end
    voxBlanks(v) = mean(GLM(v).blanks);
end
meanBetas = mean(voxBetas);
fprintf(['Saving ' hem ROI ' GLM...\n']);
% save the voxelWise mean
if ~exist('normOff','var')||normOff == 0
    eval(['save ([GLMresultsDir ''/voxelWise_GLM_' hem ROI '.mat'' ],''GLM'',''voxBetas'',''voxBlanks'')']);
else eval(['save ([GLMresultsDir ''/voxelWise_GLM_' hem ROI '_nsBetas.mat'' ],''GLM'',''voxBetas'',''voxBlanks'')']); end


