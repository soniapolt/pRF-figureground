function prfMaster_5_GLM_groupInOut_ROIacrossSubjs(expt,ROI,percCutoff)
% 10/5 edit to change voxel filtering procedure slightly
% now: first sort voxels outside of stim/too small SD, then apply
% percentage cutoff
load(['exptParams/' expt]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots ROI across subjects, rather than subject across ROI. for now, just
% voxels inside the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;
% close all;

voxPlot.ROI = ROI;
voxPlot.subjs = exptSubjs;
voxPlot.stimEdge = stimEdge;  % stim edge in ppd
voxPlot.percCutoff = percCutoff; % for each subject, restrict plots by the best X percentile of fits in each ROI

voxPlot.fitName = 'sG_multiBars';
% voxPlot.ROI ='LGN';
% voxPlot.subjs = {'160728F099','160616M090','160627F008','160721M100','160727F094','160616F007'};
voxPlot.ppd = 45.73;
voxPlot.minEccen = .251 * voxPlot.ppd; % minimum radius (to avoid edge voxels)
voxPlot.maxEccen = 4.49 * voxPlot.ppd; % maximum radius (to avoid edge voxels)
voxPlot.minSize = 0.05 * voxPlot.ppd;
% voxPlot.stimEdge = 2;  % stim edge in ppd
% voxPlot.percCutoff = 75; % for each subject, restrict plots by the best X percentile of fits in each ROI

%fMRIdir = '/users/tong_processor/Desktop/fMRI';
%expt = 'LGNedgeAll';

%plotting:
% 'aligned' 'congr-noGap' 'incongr-noGap' 'congr-gap' 'incongr-gap'
%condColors = {[140 140 140]/255; [136 188 32]/255;  [255 100 100]/255; [0 100 0]/255;[178 34 34]/255};  % grey, green,  dark green,red,dark red
lineColor = [76 0 153]/255;
barColor = [192 192 192]/255;
fontSize = 11;
titleSize = 16;
legendSize = 9;
TCfig = figure; set(gcf,'color',[1 1 1],'Units', 'Normalized', 'OuterPosition', [.1 .1 .9 .9]);
betasFig = figure; set(gcf,'color',[1 1 1],'Units', 'Normalized', 'OuterPosition', [.1 .1 .9 .9]);
hems = {'lh','rh'};



for s = 1:length(voxPlot.subjs)
    subj = voxPlot.subjs{s};
    condsIn = struct('name',condNames,'allTCs',[],'timeCourse',[],'stErrorTC',[],'offsetMeans',[]);
    condsIn = condsIn;
    voxCount = 0;
    voxIn = 0; voxOut = 0;
    
    clear bilatData; clear bilatGLM; clear bilatFits;
    
    for n = 1:length(hems)
        
        clear PRFfits; clear GLM; clear edgeData;
        fitDir = [fMRIdir '/' expt '/' subj '/matlabAnalysis/modelFit/' voxPlot.fitName '/'];
        fitFile = [hems{n} ROI '/fits_' voxPlot.fitName '_' hems{n} ROI '.mat'];
        load([fitDir fitFile]);
        
        GLMDir = [fMRIdir '/' expt '/' subj '/matlabAnalysis/GLMresults/'];
        GLMFile = ['voxelWise_GLM_' hems{n} ROI '.mat'];
        load([GLMDir GLMFile]);
        
        pscDir = [fMRIdir '/' expt '/' subj '/matlabAnalysis/meanPSC/'];
        dataFile = ['voxelWise_meanExpt_' hems{n} ROI '.mat'];
        load([pscDir dataFile]);
        
        %outputDir = [fMRIdir '/' expt  '/' subj '/matlabAnalysis/GLMinOut/'];if ~exist(outputDir) mkdir(outputDir); end
        
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
    
    stimVox = []; plotVox = [];
    for v = 1:length(bilatFits)
        if bilatFits(v).rad < voxPlot.maxEccen && bilatFits(v).rad > voxPlot.minEccen && bilatFits(v).SD > voxPlot.minSize % not on edge, not too small
            stimVox  = [stimVox  v];end
    end
    %length(inStim)
    voxPlot.R2cutoff = prctile([bilatFits(stimVox).r2],voxPlot.percCutoff);
    for v = 1:length(stimVox) if bilatFits(stimVox(v)).r2 > voxPlot.R2cutoff plotVox = [plotVox stimVox(v)]; end
    end
    %length(plotVox)
    clear PRFfits; clear GLM; clear PSCdata;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % proceed with good fits!                                                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for v = 1:length(plotVox)
        
        if isequal(bilatFits(plotVox(v)).volumeInd,bilatData(plotVox(v)).cond(1).volumeInd,bilatGLM(plotVox(v)).volumeInd)==0
            fprintf('UhOh! Fit/data/GLM might not come from the same voxels...\n');
            return
        end
        voxCount = voxCount+1;
        
        % distance from edge
        if bilatFits(plotVox(v)).rad/voxPlot.ppd < voxPlot.stimEdge
            in = 1;
            dist = voxPlot.stimEdge - (bilatFits(plotVox(v)).rad/voxPlot.ppd + bilatFits(plotVox(v)).SD/voxPlot.ppd);
        else
            in = 0;
            dist =  (bilatFits(plotVox(v)).rad/voxPlot.ppd - bilatFits(plotVox(v)).SD/voxPlot.ppd) - voxPlot.stimEdge;
        end
        
        if dist <0 || in ==1 % overlapping or in stim
            %if dist <0 && in ==0 % overlapping from outside stim
            % conds: 'aligned' 'congr-gap' 'congr-noGap' 'incongr-gap' 'incongr-noGap'
            voxIn = voxIn+1;
            for c = 1:length(condNames)
                condsIn(c).allTCs =[condsIn(c).allTCs; bilatData(plotVox(v)).cond(c).meanTC];
                condsIn(c).offsetMeans = [condsIn(c).offsetMeans ; bilatGLM(plotVox(v)).betas(c)];
            end
        end
    end % voxels
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(TCfig); % fig 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % timeCourse calculation
    totalTimeCourse = PSCparams.blockLength + PSCparams.dispPost + PSCparams.dispPre+1;
    
    for n = 1:length(condsIn)
        condsIn(n).meanTimeCourse = nanmean(condsIn(n).allTCs,1);
        condsIn(n).stErrorTimeCourse = nanstd(condsIn(n).allTCs)/sqrt(length(condsIn(n).allTCs));
    end
    
    % timeCourse plot
    subplot(2,round(length(voxPlot.subjs)/2),s);
    
    
    for n = 1:length(condsIn)
        %figure;
        errorbar([1:totalTimeCourse], 100*condsIn(n).meanTimeCourse, 100*condsIn(n).stErrorTimeCourse, 'Color',condColors{n}); hold on;
        %title(condMeans(n).name);
    end
    
    
    % mark off the block start
    line = vline((1+PSCparams.dispPre));
    set(line,'linewidth',.5,'linestyle',':','color',lineColor);
    hold on; line = vline((1+PSCparams.dispPre+PSCparams.blockLength));
    set(line,'linewidth',.5,'linestyle',':','color',lineColor);
    
    % basic properties
    set(gca,'box','off','FontSize',fontSize,'FontName','Arial','FontWeight','normal','TickDir','out');
    set(gcf,'color',[1 1 1],'Units', 'Normalized', 'OuterPosition', [.1 .1 .9 .9]);
    
    % titles and axes
    title(['Subj ' subj ': R2 > ' num2str(voxPlot.R2cutoff) ', ' num2str(voxCount) ' voxels)'],'Interpreter','none','fontsize',titleSize,'color',lineColor);
    
    xlabel('Time (TRs)','FontSize',titleSize)
    ylabel('Mean PSC','FontSize',titleSize)
    set(gca,'XTick',[1:1:totalTimeCourse],'fontSize',fontSize);
    set(gca,'XTickLabel',[-PSCparams.dispPre:1:PSCparams.blockLength+PSCparams.dispPost+1],'fontSize',fontSize*.75);
    
    % legend
    lg = legend(condsIn.name,'Location','Best');
    set(lg, 'fontsize',legendSize,'color','none','box','off');
    
    
    % basic properties
    set(gca,'box','off','FontSize',fontSize,'FontName','Arial','FontWeight','normal','TickDir','out');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(betasFig); % fig 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(2,round(length(voxPlot.subjs)/2),s);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fig 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % % trial triggered average calculation
    %figure;
    
    for n = 1:length(condsIn)
        condsIn(n).voxMean = nanmean(condsIn(n).offsetMeans);
        condsIn(n).voxSTE = nanstd(condsIn(n).offsetMeans)/sqrt(length(condsIn(n).offsetMeans));
    end
    
    hbar = bar([condsIn.voxMean]); hold on;
    errorbar([condsIn.voxMean],[condsIn.voxSTE],'linestyle', 'none');% Plot with errorbars
    set(hbar,'facecolor',barColor,'edgecolor','none');
    
    % labels
    xl = xticklabel_rotate(1:length(condsIn),70,{condsIn.name},'fontsize',fontSize+1, 'fontweight','bold');
    for n = 1:length(condsIn) set(xl(n),'color',condColors{n}); end
    ylabel('Beta Weights','FontSize',titleSize)
    title(['Subj ' subj ': R2 > ' num2str(voxPlot.R2cutoff) ', ' num2str(voxCount) ' voxels)'],'Interpreter','none','fontsize',titleSize,'color',lineColor);
    
    % basic properties
    set(gca,'box','off','FontSize',fontSize,'FontName','Arial','FontWeight','normal','TickDir','out');
    set(gcf,'color',[1 1 1],'Units', 'Normalized', 'OuterPosition', [.1 .1 .9 .9]);
end % subjs

% super title
figure(TCfig);
superTitle([expt ' Expt: Bilat ' voxPlot.ROI ' (' num2str(voxPlot.percCutoff) 'th percentile), voxels within/touching figure']);
%set(sp,'fontsize',22,'fontWeight','bold')

% super title
figure(betasFig);
superTitle([expt ' Expt: Bilat ' voxPlot.ROI ' (' num2str(voxPlot.percCutoff) 'th percentile), voxels within/touching figure']);
%set(sp,'fontsize',22,'fontWeight','bold')
