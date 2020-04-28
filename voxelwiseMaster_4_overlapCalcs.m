% edited for FT's comments, in which distance isn't calcuated but the 
clear all; close all;
expt = 'LGNedgeAll'; subjNums = 1:6; ROI = 'V1'; cutoff = 67; plotInfo = 1;

load(['exptParams/' expt]);
diffNums = diffNums(2:3,:);
diffNames = {diffNames{2} diffNames{3}};

% diffNums = diffNums(1,:);
% diffNames = {diffNames{1}};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eccenParams.fitName = 'sG_multiBars';
eccenParams.ROI = ROI;
eccenParams.subjs = exptSubjs(subjNums);%
if cutoff > 1 eccenParams.limitBy = 'perc'; else eccenParams.limitBy = 'R2'; end
eccenParams.binSize = .25;
eccenParams.edgeEccen = stimEdge;
eccenParams.plotPRFinfo = plotInfo;
eccenParams.ppd = 45.73;
eccenParams.whichSize = 'FWHM';

fontSize = 12; titleSize = 24;
sameColor = [25 25 112]/255;

hems = {'lh' 'rh' };%
voxCount = 0;
voxInfo = struct;
for s = 1:length(eccenParams.subjs)
    subj = eccenParams.subjs{s};
    
    % now we load in the data from both hemispheres, and threshold across
    % them, not within them (controls for somewhat large fit quality
    % differences across hems in some subjects
    for n = 1:length(hems)
        
        clear PRFfits; clear GLM;
        
        fitDir = [fMRIdir '/' expt '/' subj '/matlabAnalysis/modelFit/' eccenParams.fitName '/'];
        fitFile = [hems{n} eccenParams.ROI '/fits_' eccenParams.fitName '_' hems{n} eccenParams.ROI '.mat'];
        load([fitDir fitFile]);
        
        dataDir = [fMRIdir '/' expt '/' subj '/matlabAnalysis/GLMresults/'];
        dataFile = ['voxelWise_GLM_' hems{n} eccenParams.ROI '.mat'];
        load([dataDir dataFile]);
        
        if n ==1
            bilatFits = PRFfits;
            bilatGLM = GLM;
        elseif n ==2
            bilatFits = [bilatFits PRFfits];
            bilatGLM = [bilatGLM GLM];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % which voxels are we plotting? trim edge values, R2 cutoff                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(eccenParams.limitBy,'perc')==1
        [plotVox,eccenParams.R2cutoff] = analysisMaster_pickPRFs(bilatFits,cutoff,'perc');
        eccenParams.percCutoff=cutoff;
    else [plotVox,~] = analysisMaster_pickPRFs(bilatFits,cutoff,'R2');
        eccenParams.percCutoff=[]; eccenParams.R2cutoff = cutoff;
    end
    
    clear PRFfits; clear GLM;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % proceed with good fits!                                                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %voxInfo = struct;
    for v = 1:length(plotVox)
        voxCount = voxCount+1; % across subjects
        voxInfo.eccens(voxCount) = bilatFits(plotVox(v)).rad/eccenParams.ppd;
        voxInfo.Xs(voxCount) = bilatFits(plotVox(v)).XY(1)./eccenParams.ppd;
        voxInfo.Ys(voxCount) = bilatFits(plotVox(v)).XY(2)./eccenParams.ppd;
        voxInfo.SDs(voxCount) = bilatFits(plotVox(v)).SD/eccenParams.ppd;
        voxInfo.R2(voxCount) = bilatFits(plotVox(v)).r2;
        voxInfo.betas(voxCount,:) = bilatGLM(plotVox(v)).betas;
        
        % distance from edge
        if strcmp(eccenParams.whichSize,'FWHM') == 1
            voxSize = voxInfo.SDs(voxCount) * 1.1774;
        else voxSize = voxInfo.SDs(voxCount); end
        if voxInfo.eccens(voxCount) < eccenParams.edgeEccen
            voxInfo.distances(voxCount) = eccenParams.edgeEccen - (voxInfo.eccens(voxCount) + voxSize);
        else  voxInfo.distances(voxCount) =  (voxInfo.eccens(voxCount) - voxSize) - eccenParams.edgeEccen; end
        
        voxInfo.voxInd(voxCount) = bilatFits(plotVox(v)).volumeInd;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PLOT EFFECTS BY DISTANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; set(gcf,'color',[1 1 1],'Units', 'Normalized', 'OuterPosition', [.1 .1 .8 1]); 
binColors = {[0 0 1], [0.5490 0.5490 0.5490]}; % blue green grey[0.5333    0.7373    0.1255],
spDims = [2 2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
in = []; out = in;

in.vox = find([voxInfo.eccens < stimEdge]);
out.vox = setdiff([1:voxCount],in.vox);

    in.dist = -voxInfo.distances(in.vox);
    out.dist = -voxInfo.distances(out.vox); 
    xlab = 'Overlap With Boundary (dva)';


for d = 1:size(diffNums,1)
    in.diffs(d,:) = voxInfo.betas(in.vox,diffNums(d,1)) - voxInfo.betas(in.vox,diffNums(d,2));
    out.diffs(d,:) = voxInfo.betas(out.vox,diffNums(d,1)) - voxInfo.betas(out.vox,diffNums(d,2));  
end

% for now, a more complicated difference op has to be hard-coded here
if strcmp(expt,'LGNedgeAll') == 1
    in.diffs(d+1,:) = ((voxInfo.betas(in.vox,2) - voxInfo.betas(in.vox,4))+(voxInfo.betas(in.vox,3) - voxInfo.betas(in.vox,5)))/2;
    out.diffs(d+1,:) = ((voxInfo.betas(out.vox,2) - voxInfo.betas(out.vox,4))+(voxInfo.betas(out.vox,3) - voxInfo.betas(out.vox,5)))/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% voxels inside the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
in.binEdges = -1:eccenParams.binSize:max(in.dist);
[h,whichBin] = histc(in.dist,in.binEdges);

for i = 1:length(in.binEdges)
    for d = 1:length(diffNames)
    flagBinMembers = (whichBin == i);
    binMembers     = in.diffs(d,flagBinMembers);
    in.binMean(d,i)     = nanmean(binMembers);
    in.binSte(d,i) = nanstd(binMembers)/sqrt(length(binMembers));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% binned IN
subplot(spDims(1),spDims(2),1);
for d = 1:size(diffNums,1)
bhh(d) = boundedline(in.binEdges,in.binMean(d,:),in.binSte(d,:),'nan', 'gap','cmap',binColors{d},'alpha'); hold on;
[r,p] = corrcoef([in.dist in.dist],[in.diffs(d,:) in.diffs(d,:)]);
in.Rs(d) = r(2,1);in.Ps(d) = p(2,1);
end

l = legend(bhh,diffNames); set(l,'FontSize',20,'box','off');
xlabel(xlab); ylabel('Difference in Estimated Betas');
title('Voxels INSIDE the stimulus');vline(0,'k:');hline(0,'k:'); set(gca,'FontSize',20,'FontName','Helvetica Neue');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% scatter IN
subplot(spDims(1),spDims(2),spDims(2)+1);
for d = 1:size(diffNums,1)
bh(d) = scatter([in.dist],[in.diffs(d,:)],10,binColors{d},'filled'); hold on; s = lsline; hold on; 
legendIn{d} = [diffNames{d} ', R= ' num2str(in.Rs(d))]; set(s,'LineWidth',3);
end

l = legend(bh,legendIn); set(l,'FontSize',20,'box','off');
xlabel(xlab); ylabel('Difference in Estimated Betas');
title('Voxels INSIDE the stimulus');vline(0,'k:');hline(0,'k:'); set(gca,'FontSize',20,'FontName','Helvetica Neue');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% voxels outside the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
out.binEdges = -1:eccenParams.binSize:max(out.dist);
[h,whichBin] = histc(out.dist,out.binEdges);

for i = 1:length(out.binEdges)
    for d = 1:length(diffNames)
    flagBinMembers = (whichBin == i);
    binMembers     = out.diffs(d,flagBinMembers);
    out.binMean(d,i)     = nanmean(binMembers);
    out.binSte(d,i) = nanstd(binMembers)/sqrt(length(binMembers));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% binned OUT
subplot(spDims(1),spDims(2),2);
for d = 1:size(diffNums,1)
bh(d) = boundedline(out.binEdges,out.binMean(d,:),out.binSte(d,:),'nan', 'gap','cmap',binColors{d},'alpha'); hold on;
[r,p] = corrcoef([out.dist out.dist],[out.diffs(d,:) out.diffs(d,:)]);
out.Rs(d) = r(2,1);out.Ps(d) = p(2,1);
end

l = legend(bh,diffNames); set(l,'FontSize',20,'box','off');
xlabel(xlab); ylabel('Difference in Estimated Betas');
title('Voxels OUTSIDE the stimulus');vline(0,'k:');hline(0,'k:'); set(gca,'FontSize',20,'FontName','Helvetica Neue');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% scatter OUT
subplot(spDims(1),spDims(2),spDims(2)+2);
for d = 1:size(diffNums,1)
bhh(d) = scatter([out.dist],[out.diffs(d,:)],10,binColors{d},'filled'); hold on; s = lsline; hold on; 
legendOut{d} = [diffNames{d} ', R= ' num2str(out.Rs(d))]; set(s,'LineWidth',3);
end

l = legend(bhh,legendOut); set(l,'FontSize',20,'box','off');
xlabel(xlab); ylabel('Difference in Estimated Betas');
title('Voxels OUTSIDE the stimulus');vline(0,'k:');hline(0,'k:'); set(gca,'FontSize',20,'FontName','Helvetica Neue');

