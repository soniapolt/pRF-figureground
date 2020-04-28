function prfMaster_5_GLM_groupSortVox(expt,ROIs,percCutoff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multi-subj results based on pRF fits
% previous analysis yields condsIn and condsOut for each subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10/5 edit to change voxel filtering procedure slightly in lower analyses
% should have no effect here - just inherits choices in each subject's
% condsIn/condsOut data
% 8/11/17 edit takes the sortVox output rather than the old clunky output
% 
% clear all;
% close all;
saveThis = 0;
textThis = 0;

load(['exptParams/' expt]);

groupParams.fitName = 'sG_multiBars';
groupParams.ROIs = ROIs;
groupParams.percCutoff = percCutoff;
groupParams.subjs = exptSubjs;% 
groupParams.bins = {'in' 'bound' 'out' 'all'};

outputDir = [fMRIdir '/' expt '/groupResults/']; if ~exist(outputDir) mkdir(outputDir); end
 if textThis == 1 && exist('compNums')
    diaryName =  [fMRIdir 'expt_GLMttests_N'  num2str(length(exptSubjs)) '.txt']
    if exist(diaryName)>0 delete(diaryName); end
    eval(['diary ' diaryName]);
    diary ON;
    % to put at the top of our text file
    groupParams.subjs
    runTime = datestr(now)
end
   
%plotting:
% 'aligned' 'congr-noGap' 'incongr-noGap' 'congr-gap' 'incongr-gap'
lineColor = [76 0 153]/255;
barColor = [192 192 192]/255;
fontSize = 11;
titleSize = 16;
legendSize = 9;
for f = 1:length(groupParams.bins)
    eval([groupParams.bins{f} 'Fig = figure; ';]);
end
coverage = figure;

for r = 1:length(groupParams.ROIs)
conds = struct('name',condNames,'allTCs',[],'timeCourse',[],'stErrorTC',[],'subjMeans',[],'voxCount',0,'voxXYSD',[]);
group = struct('sort',groupParams.bins,'conds',conds);
counts = zeros(length(groupParams.ROIs),length(groupParams.bins));
for s= 1:length(groupParams.subjs)
    subj = groupParams.subjs{s};

    dataDir = [fMRIdir expt  '/' subj '/matlabAnalysis/GLMsorted/'];
    dataFile = [subj '_' groupParams.ROIs{r} '_perc' num2str(groupParams.percCutoff(r)) '_sortedGLM.mat'];
    load([dataDir dataFile]);
    for b = 1:length(groupParams.bins)
       counts(r,b) = counts(r,b)+ length(sortVox(b).conds(1).offsetMeans);  % number of voxels
    end
    
    for c = 1:length(group(b).conds)
        for b = 1:length(groupParams.bins)
            group(b).conds(c).allTCs =[group(b).conds(c).allTCs; sortVox(b).conds(c).meanTimeCourse];
            group(b).conds(c).subjMeans = [group(b).conds(c).subjMeans ; sortVox(b).conds(c).voxMean];
            group(b).conds(c).voxXYSD = [group(b).conds(c).voxXYSD; sortVox(b).conds(c).voxXYSD];
        end
    end
end % subjects

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% some basic ttests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('compNums')
for x = 1:length(compNames)
fprintf([groupParams.ROIs{r} ' ' num2str(percCutoff(r)) 'th percentile, ' num2str(counts(r,1)) ' total voxels inside ' num2str(stimEdge) 'dva figure\n']);
fprintf([compNames{x}  '\n']);
[~,P,~,STATS] = ttest(group(1).conds(compNums(x,1)).subjMeans,group(1).conds(compNums(x,2)).subjMeans);
fprintf(['p = ' num2str(P) ', tstat = ' num2str(STATS.tstat) ', df = ' num2str(STATS.df) '\n\n']);

fprintf([groupParams.ROIs{r} ' ' num2str(percCutoff(r)) 'th percentile, ' num2str(counts(r,2)) ' total voxels on ' num2str(stimEdge) 'dva boundary\n']);
fprintf([compNames{x}  '\n']);
[~,P,~,STATS] = ttest(group(2).conds(compNums(x,1)).subjMeans,group(2).conds(compNums(x,2)).subjMeans);
fprintf(['p = ' num2str(P) ', tstat = ' num2str(STATS.tstat) ', df = ' num2str(STATS.df) '\n\n']);

fprintf([groupParams.ROIs{r} ' ' num2str(percCutoff(r)) 'th percentile, ' num2str(counts(r,3)) ' total voxels outside ' num2str(stimEdge) 'dva figure\n']);
fprintf([compNames{x}  '\n']);
[~,P,~,STATS] = ttest(group(3).conds(compNums(x,1)).subjMeans,group(3).conds(compNums(x,2)).subjMeans);
fprintf(['p = ' num2str(P) ', tstat = ' num2str(STATS.tstat) ', df = ' num2str(STATS.df) '\n\n']);

fprintf([groupParams.ROIs{r} ' ' num2str(percCutoff(r)) 'th percentile, all ' num2str(counts(r,4)) ' voxels\n']);
fprintf([compNames{x}  '\n']);
[~,P,~,STATS] = ttest(group(4).conds(compNums(x,1)).subjMeans,group(4).conds(compNums(x,2)).subjMeans);
fprintf(['p = ' num2str(P) ', tstat = ' num2str(STATS.tstat) ', df = ' num2str(STATS.df) '\n\n']);
end
fprintf('\n');

end

% plots
for f = 1:length(sortParams.bins)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval(['figure(' sortParams.bins{f} 'Fig);']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% timeCourse calculation
PSCparams.blockLength = 8; PSCparams.dispPost = 4; PSCparams.dispPre = 4; %cheating
totalTimeCourse = PSCparams.blockLength + PSCparams.dispPost + PSCparams.dispPre+1;

for n = 1:length(group(f).conds)
    group(f).conds(n).timeCourse = nanmean(group(f).conds(n).allTCs,1);
    group(f).conds(n).stErrorTC = nanstd(group(f).conds(n).allTCs)/sqrt(length(group(f).conds(n).allTCs));
end

% timeCourse plot
subplot(2,length(groupParams.ROIs),r);

for n = 1:length(group(f).conds) 
    errorbar([1:totalTimeCourse], 100*group(f).conds(n).timeCourse, 100*group(f).conds(n).stErrorTC, 'Color',condColors{n}); hold on;
end


% mark off the block start
line = vline((1+PSCparams.dispPre)); 
set(line,'linewidth',.5,'linestyle',':','color',lineColor);
hold on; line = vline((1+PSCparams.dispPre+PSCparams.blockLength));
set(line,'linewidth',.5,'linestyle',':','color',lineColor);

% basic properties
set(gca,'box','off','FontSize',fontSize,'FontName','Arial','FontWeight','normal','TickDir','out');
set(gcf,'color',[1 1 1],'Units', 'Normalized', 'OuterPosition', [.1 .1  length(groupParams.ROIs)*.2 .8]);

% titles and axes
title(['Bilat ' groupParams.ROIs{r} ', R2=  ' num2str(groupParams.percCutoff(r)) 'th prctl: ' num2str(counts(r,f)) ' total voxels, Edge = ' sortParams.whichSize],'Interpreter','none','fontsize',titleSize,'color',lineColor);

xlabel('Time (TRs)','FontSize',titleSize)
ylabel('Mean PSC','FontSize',titleSize)
set(gca,'XTick',[1:1:totalTimeCourse],'fontSize',fontSize);
set(gca,'XTickLabel',[-PSCparams.dispPre:1:PSCparams.blockLength+PSCparams.dispPost+1],'fontSize',fontSize*.75);

% legend
lg = legend(group(f).conds.name,'Location','Best');
set(lg, 'fontsize',legendSize,'color','none','box','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% instim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % trial triggered average calculation
%figure;

for n = 1:length(group(f).conds)
    group(f).conds(n).mean = nanmean(group(f).conds(n).subjMeans);
    group(f).conds(n).STE = nanstd(group(f).conds(n).subjMeans)/sqrt(length(group(f).conds(n).subjMeans));
end

% bar plot
subplot(2,length(groupParams.ROIs),r+length(groupParams.ROIs))

hbar = bar([group(f).conds.mean]); hold on;
errorbar([group(f).conds.mean],[group(f).conds.STE],'linestyle', 'none');% Plot with errorbars
set(hbar,'facecolor',barColor,'edgecolor','none');

% labels
xl = xticklabel_rotate(1:length(group(f).conds),70,{group(f).conds.name},'fontsize',fontSize+1, 'fontweight','bold');
for n = 1:length(group(f).conds) set(xl(n),'color',condColors{n}); end
ylabel('Beta Weights','FontSize',titleSize)

% basic properties
set(gca,'box','off','FontSize',fontSize,'FontName','Arial','FontWeight','normal','TickDir','out');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% finally, a figure with all of the coverage for the paper
figure(coverage);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,length(groupParams.ROIs),r)
plotSize = sortParams.ppd*4.5;

if ~isempty([group(f).conds(1).voxXYSD])
for n = 1:size(group(f).conds(1).voxXYSD,1)
    if strcmp(sortParams.whichSize,'FWHM')==1
        rad = FWHM(group(f).conds(1).voxXYSD(n,3))/2; else rad = group(f).conds(1).voxXYSD(n,3); end
    c = plotCircle(group(f).conds(1).voxXYSD(n,1),group(f).conds(1).voxXYSD(n,2),rad,condColors{f},.5,'edge'); % center
    hold on;
    plot(group(f).conds(1).voxXYSD(n,1),-group(f).conds(1).voxXYSD(n,2),'.','Color',condColors{f}); hold on; % plot still needs Y flipped
end
end

plotCircle(0,0,sortParams.ppd*sortParams.stimEdge,'k',1,'edge');

xlim([-plotSize plotSize]); ylim([-plotSize plotSize]); axis square; hline(0,'k--'); vline(0,'k--');
title(['Group pRF Locations, ' groupParams.ROIs{r}],'FontSize',titleSize);
set(gca,'box','off','FontSize',fontSize,'FontName','Arial','FontWeight','normal','TickDir','out');
set(gcf,'color',[1 1 1],'Units', 'Normalized', 'OuterPosition', [.1 .1  length(groupParams.ROIs)*.2 .4]);

end 

if saveThis ==1
eval(['save ' outputDir 'group_' groupParams.ROIs{r} '_sortedVox.mat groupParams group']);end 



end% ROIs
diary OFF;

suptext = {'INSIDE STIM (no boundary overlap)' 'ON BOUNDARY (either in/out of stim)' 'OUTSIDE STIM (no boundary overlap)' 'ALL pRFs'};
for f = 1:length(sortParams.bins)
    eval(['figure(' sortParams.bins{f} 'Fig);']);
    % super title
    superTitle([expt ' - pRFs ' suptext{f} ' STIM, N= ' num2str(length(groupParams.subjs))]);
    %set(sp,'fontsize',22,'fontWeight','bold');
end