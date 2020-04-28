function voxelwiseMaster_4_subjRidge(expt,subjNums,ROI,cutoff,fitName,fitThis)
% ridge regression to estimate spatial profile of BOLD differences

fitName = 'cvRidge';
fitThis = 0; % will re-fit, even if this fitName already exists; otherwise, load existing fit and plot it

load(['exptParams/' expt]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what to plot?                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ridgeParams.fitName = 'sG_multiBars';
ridgeParams.ROI = ROI;
if cutoff > 1 ridgeParams.limitBy = 'perc'; else ridgeParams.limitBy = 'R2'; end
ridgeParams.edgeEccen = stimEdge;
ridgeParams.ppd = 45.73;
ridgeParams.plotppd = 4; % resolution for the plots
ridgeParams.mapSize = 4.5; % in ppd, the radius of the pRF mapping stim
ridgeParams.normWithinCond = 0; % subtract out a block average across the ROI to account for differences in arousal
ridgeParams.gaussSmooth = 3; % size of guassian smoothing filter
ridgeParams.kRange = [0:25:5000];
ridgeParams.plotSurf = 1; % either a surface plot or an imagesc plot
fontSize = 12; titleSize = 24;
sameColor = [25 25 112]/255;

hems = {'lh' 'rh' };%


ridgeDir = [fMRIdir '/' expt '/ridgeFits/'];
if ~exist(ridgeDir); mkdir(ridgeDir); end

for s = subjNums
    ridgeParams.subj = exptSubjs{s};
    voxCount = 0;
    
    ridgeFile = ['ridgeFit_' ridgeParams.subj '_' fitName '_' ROI '.mat'];
    if ~exist([ridgeDir ridgeFile]) || fitThis == 1
        
        clear voxInfo; clear ridgeFits; clear response;
        voxInfo = struct; ridgeFits = struct;
        
        % now we load in the data from both hemispheres, and threshold across
        % them, not within them (controls for somewhat large fit quality
        % differences across hems in some subjects
        for n = 1:length(hems)
            
            clear PRFfits; clear GLM;
            
            fitDir = [fMRIdir '/' expt '/' ridgeParams.subj '/matlabAnalysis/modelFit/' ridgeParams.fitName '/'];
            fitFile = [hems{n} ridgeParams.ROI '/fits_' ridgeParams.fitName '_' hems{n} ridgeParams.ROI '.mat'];
            load([fitDir fitFile]);
            
            dataDir = [fMRIdir '/' expt '/' ridgeParams.subj '/matlabAnalysis/GLMresults/'];
            dataFile = ['voxelWise_GLM_' hems{n} ridgeParams.ROI '.mat'];
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
        
        if strcmp(ridgeParams.limitBy,'perc')==1
            [plotVox,ridgeParams.R2cutoff] = analysisMaster_pickPRFs(bilatFits,cutoff,'perc');
            ridgeParams.percCutoff=cutoff;
        else [plotVox,~] = analysisMaster_pickPRFs(bilatFits,cutoff,'R2');
            ridgeParams.percCutoff=[]; ridgeParams.R2cutoff = cutoff;
        end
        
        clear PRFfits; clear GLM;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % proceed with good fits!                                                      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        diffs = [];
        
        % each subject's within-condition ROI average
        subjBlockBetas = [];
        for v=1:length(plotVox)
            subjBlockBetas(v,:,:) = reshape([bilatGLM(plotVox(v)).cond.blockBetas],length(bilatGLM(plotVox(v)).cond(1).blockBetas),length(condNames))';
        end
        subjBlockMeans{s} = squeeze(mean(subjBlockBetas));
        numBlocks(s) = length(subjBlockMeans{s});
        
        for v = 1:length(plotVox)
            voxCount = voxCount+1; % across subjects
            voxInfo.Xs(voxCount) = bilatFits(plotVox(v)).XY(1)./ridgeParams.ppd;
            voxInfo.Ys(voxCount) = bilatFits(plotVox(v)).XY(2)./ridgeParams.ppd;
            voxInfo.SDs(voxCount) = bilatFits(plotVox(v)).SD/ridgeParams.ppd;
            voxInfo.blockBetas{voxCount} = reshape([bilatGLM(plotVox(v)).cond.blockBetas],length(bilatGLM(plotVox(v)).cond(1).blockBetas),length(condNames))'; % size vox x cond x block
            voxInfo.voxInd(voxCount) = bilatFits(plotVox(v)).volumeInd;
            voxInfo.subj(voxCount) = s;
            
            if ridgeParams.normWithinCond ==1 % normalize by subtracting out ROI-wide average beta in each block
                voxInfo.blockBetas{voxCount} = voxInfo.blockBetas{voxCount}-subjBlockMeans{s}; end
            
            % mean betas (allows for recalculation after block normalization)
            voxInfo.meanBetas(voxCount,:) = mean(voxInfo.blockBetas{voxCount},2);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % regularized linear regression (ridge)                                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        GS = ridgeParams.mapSize*ridgeParams.plotppd;
        GSdim = -floor(GS):round(GS);
        [X,Y]=meshgrid(GSdim,GSdim);
        
        for v = 1:voxCount
            ridgeFits.coverage(v,:,:) = PRF(X,Y,voxInfo.Xs(v)*ridgeParams.plotppd,voxInfo.Ys(v)*ridgeParams.plotppd,voxInfo.SDs(v)*ridgeParams.plotppd);
        end
        
        for d = 1:size(diffNums,1)
            for v = 1:voxCount
                response(v) = voxInfo.meanBetas(v,diffNums(d,1)) - voxInfo.meanBetas(v,diffNums(d,2));
                weightedCoverage(v,:,:) = response(v)*squeeze(ridgeFits.coverage(v,:,:));
            end
            
            ridgeFits.response{d} = response;
            ridgeFits.weightedCoverage{d} = weightedCoverage;
            
            % regression time
            vectCoverage = reshape(ridgeFits.coverage,voxCount,size(ridgeFits.coverage,2)*size(ridgeFits.coverage,3));
            
            %kRange = logspace(0,4,500);
            tic
            [ridgeFits.betas{d},ridgeFits.chosenK(d),ridgeFits.allSSE(d,:)] = ridgeCV(response',vectCoverage,ridgeParams.kRange);
            toc
        
        
        ridgeFits.bb{d} = reshape(ridgeFits.betas{d},size(ridgeFits.coverage,2),size(ridgeFits.coverage,3));
        ridgeFits.subjVox = voxCount;
        
        save([ridgeDir ridgeFile],'ridgeParams','ridgeFits');
        end
    else % if this ridgeFit file already exists, just load it
        load([ridgeDir ridgeFile]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plots                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GS = ridgeParams.mapSize*ridgeParams.plotppd;
    GSdim = -floor(GS):round(GS);
    [X,Y]=meshgrid(GSdim,GSdim);
    ticks = [-floor(ridgeParams.mapSize):floor(ridgeParams.mapSize)];
    degMarks = GS+ridgeParams.plotppd*ticks;
    ridgeParams.plotSurf = 0;
    
    h = figure;
    set(gcf,'color',[1 1 1],'Units', 'Normalized', 'OuterPosition', [.1 .1 1 1]); set(gca,'FontSize',16);
    
    subplotDims = [length(diffNames),4]; % stimulus image,
    for d = 1:size(diffNums,1)
        % coverage
        subplot(subplotDims(1),subplotDims(2),1+(d-1)*subplotDims(2));
        imagesc(squeeze(sum(ridgeFits.coverage,1)));
        colorbar;set(gca,'XTick',degMarks,'XTickLabel',ticks,'YTick',degMarks,'YTickLabel',abs(ticks)); grid on; axis square;hline(GS,'w'); vline(GS,'w');
        title([diffNames{d} ' Sum of Coverage'],'FontSize',fontSize);hold on;
        plotCircle(GS,GS,ridgeParams.edgeEccen*ridgeParams.plotppd,[1 1 1],1,'edge'); % stim edge
        
        % coverage wieghted by responses
        subplot(subplotDims(1),subplotDims(2),2+(d-1)*subplotDims(2));
        imagesc(squeeze(sum(ridgeFits.weightedCoverage{d},1)));
        colorbar;set(gca,'XTick',degMarks,'XTickLabel',ticks,'YTick',degMarks,'YTickLabel',abs(ticks)); grid on; axis square;hline(GS,'w'); vline(GS,'w');
        title([diffNames{d} ' Sum of pRFs x Response'],'FontSize',fontSize);hold on;
        plotCircle(GS,GS,ridgeParams.edgeEccen*ridgeParams.plotppd,[1 1 1],1,'edge'); % stim edge
        
        
        
        % regression results
        subplot(subplotDims(1),subplotDims(2),3+(d-1)*subplotDims(2));
        %
        if ridgeParams.plotSurf == 1
        surf(double(squeeze(ridgeFits.bb{d}))); else
        imagesc(squeeze(ridgeFits.bb{d})); end
        colorbar;set(gca,'XTick',degMarks,'XTickLabel',ticks,'YTick',degMarks,'YTickLabel',abs(ticks)); grid on; axis square;hline(GS,'w'); vline(GS,'w');
        title([diffNames{d} ' Ridge Regression Results, MSE=' num2str(min(ridgeFits.allSSE(d,:))) ', Lambda=' num2str(ridgeFits.chosenK(d))],'FontSize',fontSize);hold on;
        if ridgeParams.plotSurf == 0 plotCircle(GS,GS,ridgeParams.edgeEccen*ridgeParams.plotppd,[1 1 1],1,'edge'); end%
        
        % smoothed regression results
        subplot(subplotDims(1),subplotDims(2),4+(d-1)*subplotDims(2));
        
        if ridgeParams.plotSurf == 1
        surf(imgaussfilt(squeeze(ridgeFits.bb{d}),ridgeParams.gaussSmooth));
        else
        imagesc(imgaussfilt(squeeze(ridgeFits.bb{d}),ridgeParams.gaussSmooth));
        end
        colorbar;set(gca,'XTick',degMarks,'XTickLabel',ticks,'YTick',degMarks,'YTickLabel',abs(ticks)); grid on; axis square;hline(GS,'w'); vline(GS,'w');
        title([diffNames{d} ' Smoothed Regression Results, GaussSize = ' num2str(ridgeParams.gaussSmooth)],'FontSize',fontSize);hold on;
        if ridgeParams.plotSurf == 0 hold on; plotCircle(GS,GS,ridgeParams.edgeEccen*ridgeParams.plotppd,[1 1 1],1,'edge'); end%%%
        
        
    end
    
    sp = suptitle([ expt ' ' (ridgeParams.subj) ' Ridge Regression pRF Reconstruction, ' ridgeParams.ROI ': ' num2str(ridgeFits.subjVox) ' voxels (' num2str(ridgeParams.percCutoff) 'th  percentile R^{2} )']);
    
    set(sp,'FontSize',titleSize);
    
    % plot the ridge parameter errors for each fit
    
%     
%     figure;
%     for d = 1:size(diffNums,1)
%         subplot(1,size(diffNums,1),d)
%         
%         plot(ridgeParams.kRange(2:end),ridgeFits.allSSE(d,2:end),'o-');
%         title(diffNames{d},'FontSize',fontSize);
%         xlabel('Value of Lamda'); ylabel('Cross-Validated MSE');
%     end
%     sp = suptitle([ expt ' Lamdba Estimation for Each Ridge Regression']);
%     set(sp,'FontSize',titleSize);set(gcf,'color',[1 1 1],'Units', 'Normalized', 'OuterPosition', [.1 .1 1 1]); set(gca,'FontSize',16);
%     
end