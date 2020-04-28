function prfMaster_4_GLM_sortVox(expt,subjNum,ROIs,percCutoff)
% 10/5 edit to change voxel filtering procedure slightly
% now: first sort voxels outside of stim/too small SD, then apply
% percentage cutoff
% 2/8/17 edit: adds another bin to separate voxels that are fully in the
% figure vs. on the border
% 8/11/17 edit: makes a nested struct because having ones named In/Out/Edge
% was ridiculous
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sorts and makes subj time courses based on the location of the pRF on the
% figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;
% close all;
saveThis = 1;

sortParams.ROIs = ROIs;
sortParams.percCutoff = percCutoff; % for each subject, restrict plots by the best X percentile of fits in each ROI
sortParams.whichSize = 'FWHM';
sortParams.bins = {'in' 'bound' 'out'};% 'all'};

load(['exptParams/' expt]);
subj = exptSubjs{subjNum};
numFuncs = subjFuncs{subjNum};

sortParams.fitName = 'sG_multiBars';
sortParams.stimEdge = stimEdge;  % stim edge in ppd
sortParams.ppd = 45.73;%

subjName = subj(end-3:end);

%plotting:
lineColor = [76 0 153]/255;
barColor = [192 192 192]/255;
fontSize = 11;
titleSize = 16;
legendSize = 9;
coverage = figure;

for f = 1:length(sortParams.bins)
    eval([sortParams.bins{f} 'Fig = figure; ';]);
end
hems = {'lh','rh'};

for r = 1:length(sortParams.ROIs)
    conds = struct('name',condNames,'allTCs',[],'offsetMeans',[],'voxInds',[],'voxXYSD',[]);
    sortVox = struct('sort',sortParams.bins,'conds',conds);
    voxCount = 0;
    
    clear bilatData; clear bilatGLM; clear bilatFits;
    
    for n = 1:length(hems)
        
        clear PRFfits; clear GLM; clear PSCdata;
        
        fitDir = [fMRIdir '/' expt '/' subj '/matlabAnalysis/modelFit/' sortParams.fitName '/'];
        fitFile = [hems{n} sortParams.ROIs{r} '/fits_' sortParams.fitName '_' hems{n} sortParams.ROIs{r} '.mat'];
        load([fitDir fitFile]);
        
        GLMDir = [fMRIdir '/' expt '/' subj '/matlabAnalysis/GLMresults/'];
        GLMFile = ['voxelWise_GLM_' hems{n} sortParams.ROIs{r} '.mat'];
        load([GLMDir GLMFile]);
        
        pscDir = [fMRIdir '/' expt '/' subj '/matlabAnalysis/meanPSC/'];
        dataFile = ['voxelWise_meanExpt_' hems{n} sortParams.ROIs{r} '.mat'];
        load([pscDir dataFile]);
        
        outputDir = [fMRIdir '/' expt  '/' subj '/matlabAnalysis/GLMsorted/'];if ~exist(outputDir) mkdir(outputDir); end
        % at some point, delet GLMinOut, the old dirs
        
        if n ==1
            bilatFits = PRFfits;
            bilatGLM = GLM;
            bilatData = PSCdata;
        elseif n ==2
            bilatFits = [bilatFits PRFfits];
            bilatGLM = [bilatGLM GLM];
            bilatData = [bilatData PSCdata];
        end
        
    end % hems
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % which voxels are we plotting? trim edge values, R2 cutoff                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [plotVox,sortParams.R2cutoffs(r)] = analysisMaster_pickPRFs(bilatFits,percCutoff(r),'perc');
    
    clear PRFfits; clear GLM; clear PSCdata;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % proceed with good fits!                                                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for v = 1:length(plotVox)
        
        if isequal(bilatFits(plotVox(v)).volumeInd,bilatData(plotVox(v)).cond(1).volumeInd,bilatGLM(plotVox(v)).volumeInd)==0
            bilatFits(plotVox(v)).volumeInd
            bilatData(plotVox(v)).cond(1).volumeInd
            bilatGLM(plotVox(v)).volumeInd
            fprintf('UhOh! Fit/data/GLM might not come from the same voxels...\n');
            return
        end
        
        voxCount = voxCount+1;
        
        if strcmp(sortParams.whichSize,'FWHM') == 1
            pRFedge =  0.5*FWHM(bilatFits(plotVox(v)).SD)/sortParams.ppd;
        else pRFedge =  bilatFits(plotVox(v)).SD/sortParams.ppd; end
        
        
        % distance from edge
        if bilatFits(plotVox(v)).rad/sortParams.ppd < sortParams.stimEdge
            in = 1;
            dist = sortParams.stimEdge - (bilatFits(plotVox(v)).rad/sortParams.ppd + pRFedge);
            
        else
            in = 0;
            dist =  (bilatFits(plotVox(v)).rad/sortParams.ppd - pRFedge) - sortParams.stimEdge;
        end
        
        % sortParams.bins = {'in' 'bound' 'out' 'all'};
        if dist >0 && in == 1 % in stim, not on boundary
            sorted = 1;
        elseif dist > 0 && in == 0 % out of stim
            sorted = 3;
        elseif dist <0  % on the boundary
            sorted = 2;
        end
        
        % conds: 'aligned' 'congr-gap' 'congr-noGap' 'incongr-gap' 'incongr-noGap'
        for c = 1:length(condNames)
            for s = [sorted length(sortParams.bins)] % always include a bin with all the data for convenience
                sortVox(s).conds(c).allTCs =[sortVox(s).conds(c).allTCs; bilatData(plotVox(v)).cond(c).meanTC];
                sortVox(s).conds(c).offsetMeans = [sortVox(s).conds(c).offsetMeans ; bilatGLM(plotVox(v)).betas(c)];
                sortVox(s).conds(c).voxInds = [sortVox(s).conds(c).voxInds; plotVox(v)];
                sortVox(s).conds(c).voxXYSD = [sortVox(s).conds(c).voxXYSD; [bilatFits(plotVox(v)).XY bilatFits(plotVox(v)).SD]];%
            end
        end
    end % voxels
    
    for f = 1:length(sortParams.bins)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        eval(['figure(' sortParams.bins{f} 'Fig);']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % timeCourse calculation
        totalTimeCourse = PSCparams.blockLength + PSCparams.dispPost + PSCparams.dispPre+1;
        
        for n = 1:length(sortVox(f).conds)
            if size(sortVox(f).conds(n).allTCs,1)>0
                sortVox(f).conds(n).meanTimeCourse = nanmean(sortVox(f).conds(n).allTCs,1);
                sortVox(f).conds(n).stErrorTimeCourse = nanstd(sortVox(f).conds(n).allTCs,1)/sqrt(length(sortVox(f).conds(n).allTCs));
            else sortVox(f).conds(n).meanTimeCourse  = zeros(1,totalTimeCourse); end
            if size(sortVox(f).conds(n).allTCs,1)<2
                sortVox(f).conds(n).stErrorTimeCourse = zeros(1,totalTimeCourse); end
        end
        
        % timeCourse plot
        subplot(2,length(sortParams.ROIs),r);
        % basic properties
        
        
        for n = 1:length(conds)
            errorbar([1:totalTimeCourse], 100*sortVox(f).conds(n).meanTimeCourse, 100*sortVox(f).conds(n).stErrorTimeCourse, 'Color',condColors{n}); hold on;
        end
        
        % mark off the block start
        line = vline((1+PSCparams.dispPre));
        set(line,'linewidth',.5,'linestyle',':','color',lineColor);
        hold on; line = vline((1+PSCparams.dispPre+PSCparams.blockLength));
        set(line,'linewidth',.5,'linestyle',':','color',lineColor);
        
        set(gca,'box','off','FontSize',fontSize,'FontName','Arial','FontWeight','normal','TickDir','out');
        set(gcf,'color',[1 1 1],'Units', 'Normalized', 'OuterPosition', [.1 .1  length(sortParams.ROIs)*.2 .8]);
        % titles and axes
        title(['Bilat ' sortParams.ROIs{r} ' (' num2str(length(sortVox(f).conds(1).voxInds)) ' voxels, R2 > ' num2str(sortParams.R2cutoffs(r)) ', ' num2str(sortParams.percCutoff(r)) 'th perc)'],'Interpreter','none','fontsize',titleSize,'color',lineColor);
        
        xlabel('Time (TRs)','FontSize',titleSize)
        ylabel('Mean PSC','FontSize',titleSize)
        set(gca,'XTick',[1:1:totalTimeCourse],'fontSize',fontSize);
        set(gca,'XTickLabel',[-PSCparams.dispPre:1:PSCparams.blockLength+PSCparams.dispPost+1],'fontSize',fontSize*.75);
        
        
        % legend
        lg = legend(sortVox(f).conds.name,'Location','Best');
        set(lg, 'fontsize',legendSize,'color','none','box','off');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % trial triggered average calculation
        
        for n = 1:length(sortVox(f).conds)
            sortVox(f).conds(n).voxMean = nanmean(sortVox(f).conds(n).offsetMeans);
            if length(sortVox(f).conds(n).offsetMeans)>1
                sortVox(f).conds(n).voxSTE = nanstd(sortVox(f).conds(n).offsetMeans)/sqrt(length(sortVox(f).conds(n).offsetMeans));
            else sortVox(f).conds(n).voxSTE = 0; end
        end
        
        
        % bar plot
        subplot(2,length(sortParams.ROIs),r+length(sortParams.ROIs))
        
        hbar = bar([sortVox(f).conds.voxMean]); hold on;
        errorbar([sortVox(f).conds.voxMean],[sortVox(f).conds.voxSTE],'linestyle', 'none');% Plot with errorbars
        set(hbar,'facecolor',barColor,'edgecolor','none');
        
        
        % labels
        xl = xticklabel_rotate(1:length(sortVox(f).conds),70,{sortVox(f).conds.name},'fontsize',fontSize+1, 'fontweight','bold');
        for n = 1:length(sortVox(f).conds) set(xl(n),'color',condColors{n}); end
        ylabel('Beta Weights','FontSize',titleSize)
        
        % basic properties
        set(gca,'box','off','FontSize',fontSize,'FontName','Arial','FontWeight','normal','TickDir','out');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% finally, a figure with all of the coverage for the paper
        figure(coverage);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(1,length(sortParams.ROIs),r)
        plotSize = sortParams.ppd*4.5;
        
        if ~isempty([sortVox(f).conds(1).voxInds])
            for v = 1:length([sortVox(f).conds(1).voxInds])
                n = sortVox(f).conds(1).voxInds(v);
                if strcmp(sortParams.whichSize,'FWHM')==1
                    rad = FWHM(bilatFits(n).SD)/2; else rad = bilatFits(n).SD; end
                hold on;
                c = plotCircle(bilatFits(n).XY(1),bilatFits(n).XY(2),rad,condColors{f},.5,'edge'); % center
                hold on;
                plot(bilatFits(n).XY(1),-bilatFits(n).XY(2),'.','Color',condColors{f}); hold on; % plot still needs Y flipped
            end
        end
        
        xlim([-plotSize plotSize]); ylim([-plotSize plotSize]); axis square; hline(0,'k--'); vline(0,'k--');
        title(['Individual pRF Locations, ' subjName ' ' sortParams.ROIs{r}],'FontSize',titleSize);
        set(gca,'box','off','FontSize',fontSize,'FontName','Arial','FontWeight','normal','TickDir','out');
        set(gcf,'color',[1 1 1],'Units', 'Normalized', 'OuterPosition', [.1 .1  length(sortParams.ROIs)*.2 .4]);
        
        plotCircle(0,0,sortParams.ppd*sortParams.stimEdge,'k',1,'edge');
    end
    
    if saveThis ==1
        eval(['save ([outputDir subj ''_'' sortParams.ROIs{r} ''_perc'' num2str(sortParams.percCutoff(r)) ''_sortedGLM.mat'' ],''sortVox'',''sortParams'');']); end
end % ROI


suptext = {'INSIDE STIM (no boundary overlap)' 'ON BOUNDARY (either in/out of stim)' 'OUTSIDE STIM (no boundary overlap)' 'ALL VOX'};
for f = 1:length(sortParams.bins)
    eval(['figure(' sortParams.bins{f} 'Fig);']);
    % super title
    superTitle([expt ' - pRFs ' suptext{f} ', Subj. ' subjName ', pRF edge = ' sortParams.whichSize]);
    %set(sp,'fontsize',22,'fontWeight','bold');
end
