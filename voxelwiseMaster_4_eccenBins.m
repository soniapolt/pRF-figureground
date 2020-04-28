%function voxelwiseMaster_4_eccenBins(expt,subjNums,ROI,cutoff,plotInfo)
% prf bin properties and binned plots

clear all; close all;
expt = 'LGNedgeAll'; subjNums = 1:6; ROI = 'V1'; cutoff = 67; plotInfo = 1;
%expt = 'LGNedgeSmall'; subjNums = 1; ROI = 'V1'; cutoff = 50; plotInfo = 1;
%expt = 'LGNedge3Pilot'; subjNums = 1:3; ROI = 'LGN'; cutoff = 50; plotInfo = 1;

load(['exptParams/' expt]);
fMRIdir = dirOf(pwd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eccenParams.fitName = 'sG_multiBars';
eccenParams.ROI = ROI;
eccenParams.subjs = exptSubjs(subjNums);%
if cutoff > 1 eccenParams.limitBy = 'perc'; else eccenParams.limitBy = 'R2'; end
eccenParams.binSize = .25;
eccenParams.edgeEccen = stimEdge;
eccenParams.plotPRFinfo = 0;
eccenParams.whichSize = 'FWHM'; % 'FWHM' or 'SD'
eccenParams.ppd = 45.73;
eccenParams.minVox = 10; % don't plot bins with fewer than X voxels

fontSize = 12; titleSize = 24;
sameColor = [25 25 112]/255;

hems = {'lh' 'rh' };%
voxCount = 0;
voxInfo = struct;
subjVox = []; subCuts = [];
for s = 1:length(eccenParams.subjs)
    subj = eccenParams.subjs{s}
    
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
    subjCuts(s) = eccenParams.R2cutoff;
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
        voxInfo.FWHMs(voxCount) = bilatFits(plotVox(v)).SD/eccenParams.ppd * 1.1774;
        voxInfo.R2(voxCount) = bilatFits(plotVox(v)).r2;
        voxInfo.betas(voxCount,:) = bilatGLM(plotVox(v)).betas;
        eval(['voxInfo.size(voxCount) = voxInfo.' eccenParams.whichSize 's(voxCount);']);
        
        % distance from edge
        if voxInfo.eccens(voxCount) < eccenParams.edgeEccen
            voxInfo.distances(voxCount) = eccenParams.edgeEccen - (voxInfo.eccens(voxCount) + voxInfo.SDs(voxCount));
        else  voxInfo.distances(voxCount) =  (voxInfo.eccens(voxCount) - voxInfo.SDs(voxCount)) - eccenParams.edgeEccen; end
        
        voxInfo.voxInd(voxCount) = bilatFits(plotVox(v)).volumeInd;
    end
    subjVox(s) = length(plotVox);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chop up data by bins                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%voxCount = length(plotVox);

% binning calc
eccenParams.bins = [0:eccenParams.binSize:4.5];

% to check
a = struct('name',condNames,'vox',[]);b = struct('name',diffNames,'vals',[]);
bins(1:length(eccenParams.bins)-1) = struct('cond',a,'diff',b,'size',[],'r2',[],'dist',[],'distProp',[]);
for n= 1:length(voxInfo.eccens) % voxel
    for m = 2:length(eccenParams.bins) % bin
        if voxInfo.eccens(n) >= eccenParams.bins(m-1) && voxInfo.eccens(n) < eccenParams.bins(m)
            
            for c = 1:length(condNames)
                bins(m-1).cond(c).vox = [bins(m-1).cond(c).vox  voxInfo.betas(n,c)];
            end
            
            bins(m-1).size = [bins(m-1).size voxInfo.size(n)];
            bins(m-1).r2 = [bins(m-1).r2 voxInfo.R2(n)];
            
            bins(m-1).dist = [bins(m-1).dist voxInfo.distances(n)];
            if voxInfo.distances(n) <= 0
                bins(m-1).distProp = [bins(m-1).distProp 1];
            else bins(m-1).distProp = [bins(m-1).distProp 0];end
            
        end
    end
end

for n = 1:length(bins)
    bins(n).center = eccenParams.bins(n+1)-eccenParams.binSize/2;
    
    % prf properties
    bins(n).sizeMean = nanmean(bins(n).size);
    bins(n).sizeSte = nanstd(bins(n).size)/sqrt(length(bins(n).size));
    bins(n).numVox = length(bins(n).size);
    bins(n).r2mean = nanmean(bins(n).r2);
    bins(n).r2ste = nanstd(bins(n).r2)/sqrt(length(bins(n).r2));
    
    bins(n).distMean = nanmean(bins(n).dist);
    bins(n).distSte = nanstd(bins(n).dist)/sqrt(length(bins(n).dist));
    bins(n).distProp = mean(bins(n).distProp);

    % summary of each condition
    for c = 1:length(condNames)
        bins(n).cond(c).mean = nanmean(bins(n).cond(c).vox);
        bins(n).cond(c).ste = nanstd(bins(n).cond(c).vox)./sqrt(length(bins(n).cond(c).vox));
    end
    
    % difference calculations
    for d = 1:length(diffNums)
        bins(n).diff(d).vals = bins(n).cond(diffNums(d,1)).vox-bins(n).cond(diffNums(d,2)).vox;
        bins(n).diff(d).mean = nanmean(bins(n).diff(d).vals);
        bins(n).diff(d).ste = nanstd(bins(n).diff(d).vals)/sqrt(length(bins(n).diff(d).vals));
    end
    
    % for now, a more complicated difference op has to be hard-coded here
    if strcmp(expt,'LGNedgeAll') == 1
        diffNames{d+1} = 'Gap Effect';
        bins(n).diff(d+1).vals = ((bins(n).cond(2).vox-bins(n).cond(4).vox)+(bins(n).cond(3).vox-bins(n).cond(5).vox))/2;
        bins(n).diff(d+1).mean = nanmean(bins(n).diff(d+1).vals);
        bins(n).diff(d+1).ste = nanstd(bins(n).diff(d+1).vals)/sqrt(length(bins(n).diff(d+1).vals));
    end
    
    % NaN out bins with fewer than X voxels
    if bins(n).numVox < eccenParams.minVox
     [bins(n).diff.mean] = deal(NaN);
     [bins(n).diff.ste] = deal(NaN);
     [bins(n).cond.mean] = deal(NaN);
     [bins(n).cond.ste] = deal(NaN); 
    end
    
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% plots o'clock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numPlots = 1+length(diffNames);
if numPlots>4 numRows = 2; else numRows = 1; end
numCols = round(numPlots/numRows);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DIFFERENCE PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; set(gcf,'color',[1 1 1],'Units', 'Normalized', 'OuterPosition', [.1 .1 .6 1]); set(gca,'FontSize',16);
set(gcf,'renderer','Painters');

for d = 1:length(diffNames)
    subplot(numRows,numCols,1+d);
    
    plotDiffs = []; plotErrs = [];
    for b = 1:length(bins)
        plotDiffs = [plotDiffs bins(b).diff(d).mean];
        plotErrs = [plotErrs bins(b).diff(d).ste]; end
    hbar = bar([bins.center], [plotDiffs]); hold on;
    eb = errorbar([bins.center], [plotDiffs], [plotErrs],'linestyle', 'none');% Plot with errorbars
    errorbar_fix(eb);
    set(hbar,'facecolor',condColors{d},'edgecolor','none');
    title(diffNames{d});
    xlabel('Binned Voxel Eccentricity (dva)'); ylabel('Difference in (Beta Weight)');
    axis square;
    vline(stimEdge,'k:');
    set(gca,'box','off','FontSize',fontSize,'FontName','Arial','FontWeight','normal','TickDir','out');
end

if strcmp(eccenParams.limitBy,'perc')==1
    suptitle([eccenParams.ROI ': ' num2str(voxCount) ' voxels (each subject''s ' num2str(eccenParams.percCutoff) 'th percentile)']);
else
    suptitle([eccenParams.ROI ': ' num2str(voxCount) ' voxels (R2 > ' num2str(eccenParams.R2cutoff) ')']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CONDITIONS LINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(numRows,numCols,1);

for c = 1:3%length(condNames)
    means = []; errors = [];
    for b = 1:length(bins)
        means = [means bins(b).cond(c).mean];
        errors = [errors bins(b).cond(c).ste];
    end
    %errorbar([bins.center], means,  errors,'Color',condColors{c}); hold on;
    shadedErrorBar([bins.center], means, errors ,condColors{c},'lineprops',{'linewidth',2,'color',condColors{c}});
end

xlabel('Max Binned Eccentricity (dva)'); ylabel('Beta Weight');
xlabel('Max Binned Eccentricity (dva)'); ylabel('Beta Weight');title('No Gap');
%legend(condNames,'fontsize',fontSize,'location','best','color','none','box','off');
axis square;vline(stimEdge,'k:');
set(gca,'box','off','FontSize',fontSize,'FontName','Arial','FontWeight','normal','TickDir','out');

if strcmp(eccenParams.limitBy,'perc')==0
    if length(eccenParams.subjs)>1
        suptitle([eccenParams.ROI ': ' num2str(voxCount) ' voxels across N = ' num2str(length(eccenParams.subjs)) ' (R2 > ' num2str(eccenParams.R2cutoff) '), Bin Size = ' num2str(eccenParams.binSize) ' dva']);
    else
        suptitle([eccenParams.ROI ': ' num2str(voxCount) ' voxels, Subj = ' eccenParams.subjs{1} ' (R2 > ' num2str(eccenParams.R2cutoff) '), Bin Size = ' num2str(eccenParams.binSize) ' dva']);
    end
else
    if length(eccenParams.subjs)>1
        suptitle([eccenParams.ROI ': ' num2str(voxCount) ' voxels across N = ' num2str(length(eccenParams.subjs)) ' (' num2str(eccenParams.percCutoff) 'th percentile), Bin Size = ' num2str(eccenParams.binSize) ' dva']);
    else
        suptitle([eccenParams.ROI ': ' num2str(voxCount) ' voxels, Subj = ' eccenParams.subjs{1} ' (' num2str(eccenParams.percCutoff) 'th percentile), Bin Size = ' num2str(eccenParams.binSize) ' dva']);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PRF PROPERTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if eccenParams.plotPRFinfo ==1
    figure;set(gcf,'color',[1 1 1],'Units', 'Normalized', 'OuterPosition', [.1 .1 .9 .9]); set(gca,'FontSize',16);
    
    % size by eccen
    subplot(2,3,1);
    
    scatter([voxInfo.eccens],voxInfo.FWHMs,30);l = lsline; set(l,'LineWidth',3); hold on;
    xlabel('Radius (dva)'); ylabel('Center FWHM (dva)'); title('Center FWHM by Eccentricity');ylim([0 5]);
    axis square;
    
    xlabel('Max Binned Eccentricity (dva)'); ylabel('Width of PRF (FWHM in dva)');
    title('pRF Sizes Accross Eccen');
    
    % voxels per bin
    subplot(2,3,2);
    hbar = bar([bins.center],[bins.numVox]); hold on;
    xlabel('Max Binned Eccentricity (dva)'); ylabel('Count of Voxels');
    title('Number of Voxels in Each Bin');
    
    % R2 in each bin
    subplot(2,3,3);
    hbar = bar([bins.center],[bins.r2mean]); hold on;
    errorbar([bins.center],[bins.r2mean],[bins.r2ste],'linestyle', 'none');% Plot with errorbars
    xlabel('Max Binned Eccentricity (dva)'); ylabel('Mean R2'); ylim([min(eccenParams.R2cutoff(:))-.05 1]);
    title('Mean R2 in Each Bin');
    
    % coverage across all voxels plotted
    subplot(2,3,4);
    for v = 1:voxCount
        f = plotCircle(voxInfo.Xs(v),voxInfo.Ys(v),1.774*voxInfo.SDs(v)/2,[0 0 1],.5,'edge'); % FWHM
        set(f,'Linewidth',1); hold on;
        plot(voxInfo.Xs(v),-voxInfo.Ys(v),'k.');hold on;
    end
    % add outer edge & stim location
    f = plotCircle(0,0,4.5,[0 0 0],1,'edge'); hold on; set(f,'LineStyle','--','Linewidth',1.5);
    f = plotCircle(0,0,eccenParams.edgeEccen,[1 0 0],1,'edge'); hold on; set(f,'Linewidth',1.5);
    plotSize = 4.5;
    xlim([-plotSize plotSize]); ylim([-plotSize plotSize]); axis square; hline(0,'k--'); vline(0,'k--');
    title(['Coverage Across All ' (length(voxCount)) ' Voxels (FWHM)']);
    
    
    % average distance to the edge in each bin
    subplot(2,3,5);
    errorbar([bins.center], [bins.distMean], [bins.distSte]);
    xlabel('Max Binned Eccentricity (dva)'); ylabel('Distance to Edge');
    title('pRF Distance-to-Stim-Edge Accross Eccen');
    
    % percent of voxels that overlap
    subplot(2,3,6);
    hbar = bar([bins.center],[bins.distProp]); hold on;
    xlabel('Max Binned Eccentricity (dva)'); ylabel('Proportion of Voxels');
    title(['Proportion of Voxels that Overlap Stim at 1' eccenParams.whichSize]);
    
    if length(eccenParams.subjs)>1
        suptitle([eccenParams.ROI ': ' num2str(voxCount) ' voxels across N = ' num2str(length(eccenParams.subjs)) ' (' num2str(eccenParams.percCutoff) 'th percentile), Bin Size = ' num2str(eccenParams.binSize) ' dva']);
    else
        suptitle([eccenParams.ROI ': ' num2str(voxCount) ' voxels, Subj = ' eccenParams.subjs{1} ' (' num2str(eccenParams.percCutoff) 'th percentile), Bin Size = ' num2str(eccenParams.binSize) ' dva']);
    end
end

subjVox
mean(subjVox)
std(subjVox)

subjCuts
mean(subjCuts)