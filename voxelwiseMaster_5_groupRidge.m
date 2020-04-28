% loads in individual-subject ridge regressions, zscores the images,
% averages them.

clear all; %close all;
expt = 'LGNedgeAll'; subjNums = 1:6; ROI = 'V3';

fitName = 'cvRidge_67';

load(['exptParams/' expt]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupRidge.ROI = ROI;
groupRidge.zNorm = 1; % z-normalize each subject's reconstructed image
groupRidge.gaussSmooth = 3; % size of guassian smoothing filter
fontSize = 12; titleSize = 24;
sameColor = [25 25 112]/255;
groupRidge.edgeEccen = stimEdge;
groupRidge.ppd = 45.73;
groupRidge.plotppd = 4; % resolution for the plots
groupRidge.mapSize = 4.5; % in ppd, the radius of the pRF mapping stim


h = figure;
set(gcf,'color',[1 1 1],'Units', 'Normalized', 'OuterPosition', [.1 .1 1 1]); set(gca,'FontSize',16);
subplotDims = [length(diffNames),4];

for d = 1:size(diffNums,1)
    groupRidge.coverage = []; groupRidge.weightedCoverage = groupRidge.coverage; groupRidge.allBBs = groupRidge.coverage;
    for s = subjNums
        subj = exptSubjs{s};
        
        ridgeDir = [fMRIdir '/' expt '/ridgeFits/'];
        ridgeFile = ['ridgeFit_' subj '_' fitName '_' groupRidge.ROI '.mat'];
        if ~exist([ridgeDir ridgeFile])
            fprintf(['Uh oh, subject ' subj ' does not have individual reconstruction!\n']);
            break
        else
            load([ridgeDir ridgeFile]);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % aggreggate coverage & fits                                                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         if size(ridgeFits.betas)==1
        %             bb = reshape(ridgeFits.betas{1},size(ridgeFits.coverage,2),size(ridgeFits.coverage,3));
        %             groupRidge.weightedCoverage(s,:,:)  = squeeze(mean(ridgeFits.weightedCoverage{1}));
        %         else
        %             bb = reshape(ridgeFits.betas,size(ridgeFits.coverage,2),size(ridgeFits.coverage,3));
        %             groupRidge.weightedCoverage(s,:,:)  = squeeze(mean(ridgeFits.weightedCoverage));
        %         end%
        
        bb = ridgeFits.bb{d};
        groupRidge.coverage(s,:,:) = squeeze(mean(ridgeFits.coverage));
        groupRidge.weightedCoverage(s,:,:)  = squeeze(mean(ridgeFits.weightedCoverage{d}));
        
        
        
        if groupRidge.zNorm == 1
            bb = (bb-mean(bb(:)))/std(bb(:));
        end
        
        groupRidge.allBBs(s,:,:)  = bb;
    end
    
    groupRidge.bb = squeeze(mean(groupRidge.allBBs));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plots - group                                                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GS = ridgeParams.mapSize*ridgeParams.plotppd;
    ticks = [-floor(groupRidge.mapSize):floor(groupRidge.mapSize)];
    degMarks = GS+groupRidge.plotppd*ticks;
    
    % coverage
    subplot(subplotDims(1),subplotDims(2),1+(d-1)*subplotDims(2));
    imagesc(squeeze(sum(groupRidge.coverage,1)));
    colorbar;set(gca,'XTick',degMarks,'XTickLabel',ticks,'YTick',degMarks,'YTickLabel',abs(ticks)); grid on; axis square;hline(GS,'w'); vline(GS,'w');
    title([diffNames{d} ' Sum of Coverage'],'FontSize',fontSize);hold on;
    plotCircle(GS,-GS,groupRidge.edgeEccen*groupRidge.plotppd,[1 1 1],1,'edge'); % stim edge
    if strcmp(diffNames{d},'Figure (Gap)') == 1
        plotCircle(GS,-GS,(groupRidge.edgeEccen+.5)*groupRidge.plotppd,[1 1 1],1,'edge'); % stim edge
    end
    % coverage wieghted by responses
    subplot(subplotDims(1),subplotDims(2),2+(d-1)*subplotDims(2));
    imagesc(squeeze(sum(groupRidge.weightedCoverage,1)));
    colorbar;set(gca,'XTick',degMarks,'XTickLabel',ticks,'YTick',degMarks,'YTickLabel',abs(ticks)); grid on; axis square;hline(GS,'w'); vline(GS,'w');
    title([diffNames{d} ' Sum of pRFs x Response'],'FontSize',fontSize);hold on;
    plotCircle(GS,-GS,groupRidge.edgeEccen*groupRidge.plotppd,[1 1 1],1,'edge'); % stim edge
    if strcmp(diffNames{d},'Figure (Gap)') == 1
        plotCircle(GS,-GS,(groupRidge.edgeEccen+.5)*groupRidge.plotppd,[1 1 1],1,'edge'); % stim edge
    end
    % regression results
    subplot(subplotDims(1),subplotDims(2),3+(d-1)*subplotDims(2));
    imagesc(groupRidge.bb);
    colorbar;set(gca,'XTick',degMarks,'XTickLabel',ticks,'YTick',degMarks,'YTickLabel',abs(ticks)); grid on; axis square;hline(GS,'w'); vline(GS,'w');
    title([diffNames{d} ' Ridge Regression Average Results, Z-Norm = ' num2str(groupRidge.zNorm)],'FontSize',fontSize);hold on;
    plotCircle(GS,-GS,groupRidge.edgeEccen*groupRidge.plotppd,[1 1 1],1,'edge'); % stim edge
    if strcmp(diffNames{d},'Figure (Gap)') == 1
        plotCircle(GS,-GS,(groupRidge.edgeEccen+.5)*groupRidge.plotppd,[1 1 1],1,'edge'); % stim edge
    end
    % smoothed regression results
    subplot(subplotDims(1),subplotDims(2),4+(d-1)*subplotDims(2));
    imagesc(imgaussfilt(groupRidge.bb,groupRidge.gaussSmooth));
    colorbar;set(gca,'XTick',degMarks,'XTickLabel',ticks,'YTick',degMarks,'YTickLabel',abs(ticks)); grid on; axis square;hline(GS,'w'); vline(GS,'w');
    title([diffNames{d} ' Smoothed Regression Results, GaussSize = ' num2str(groupRidge.gaussSmooth)],'FontSize',fontSize);hold on;
    plotCircle(GS,-GS,groupRidge.edgeEccen*groupRidge.plotppd,[1 1 1],1,'edge'); % stim edge
    if strcmp(diffNames{d},'Figure (Gap)') == 1
        plotCircle(GS,-GS,(groupRidge.edgeEccen+.5)*groupRidge.plotppd,[1 1 1],1,'edge'); % stim edge
    end
    
end
sp = suptitle([ expt ' Ridge Regression pRF Reconstruction, ' groupRidge.ROI ', N=' num2str(length(subjNums))]);
set(sp,'FontSize',titleSize);
