function [plotVox,R2cutoff] = analysisMaster_pickPRFs(bilatFits,cutoff,limitBy)
% unified way to pick pRFs from existing fits & trimming
% cutoff can either be a percentile or an r2 value, indicated by limitBy
% 'perc' or 'R2'

trim.ppd = 45.73; % value from display, matched to fitting code
trim.minEccen = .251 * trim.ppd; % minimum radius (to avoid edge voxels)
trim.maxEccen = 4.49 * trim.ppd; % maximum radius (to avoid edge voxels)
trim.minSize = 0.05 * trim.ppd; % minimum size (to avoid zero-response estimates)

stimVox = []; plotVox = [];
    for v = 1:length(bilatFits)
            if bilatFits(v).rad < trim.maxEccen && bilatFits(v).rad > trim.minEccen && bilatFits(v).SD > trim.minSize % not on edge, not too small
             stimVox  = [stimVox  v];end
    end
    if strcmp(limitBy,'perc')==1
    R2cutoff = prctile([bilatFits(stimVox).r2],cutoff);
    else R2cutoff = cutoff; end
    for v = 1:length(stimVox) if bilatFits(stimVox(v)).r2 > R2cutoff plotVox = [plotVox stimVox(v)]; end 
    end

end

