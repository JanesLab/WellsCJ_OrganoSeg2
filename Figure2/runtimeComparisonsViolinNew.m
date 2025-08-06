x49 = cell(1,98);
clf
filenames = {dir().name};
for filename = filenames
    if contains(filename, "final_runtime")
        
        load(filename{1})

    end
end

for i = 1:98
    if i <= 49
        x49{i} = 'OrganoSeg';
    else
        x49{i} = 'OrganoSeg2';
    end
end
x49 = categorical(x49)';

ySeg = [timesSegment'; timesSegmentNew'];
yexpDef = [timesExport'; timesExportNew'];
yexpAll = [timesExportAll'; timesExportNewAll'];
yImage = [timesChange'; timesChangeNew'];





%% 


close all
figure(1)
h = violinplot(log10(ySeg), x49,'ShowBox', false, 'ViolinColor',  ([0,0,200]/255));
parentAxes = h.Parent;
stylegraph(parentAxes)
parentAxes.YLim = [-1.5 1.5];
[pSeg, h, statsSeg] = ranksum(timesSegment, timesSegmentNew);
title('Segmentation')
ylabel("Log Runtime (s)")


figure(2)
h = violinplot(log10(yexpDef), x49,'ShowBox', false, 'ViolinColor', ([0,0,200]/255));
parentAxes = h.Parent;

stylegraph(parentAxes)
parentAxes.YLim = [-1.5 1.5];
[pExpDef, h, statsExpDef] = ranksum(timesExport, timesExportNew);
title('Export Default Metrics')

ylabel("Log Runtime (s)")


figure(3)
h = violinplot(log10(yexpAll), x49,'ShowBox', false, 'ViolinColor', ([0,0,200]/255));
parentAxes = h.Parent;
stylegraph(parentAxes)
parentAxes.YLim = [-1.5 1.5];
[pExpAll, h, statsExpAll] = ranksum(timesExportAll, timesExportNewAll);
title('Export All Metrics')
ylabel("Log Runtime (s)")


figure(4)
h = violinplot(log10(yImage), x49,'ShowBox', false, 'ViolinColor', ([0,0,200]/255));
parentAxes = h.Parent;
stylegraph(parentAxes)
parentAxes.YLim = [-1.5 1.5];
[pImage, h, statsImage] = ranksum(timesChange, timesChangeNew);
title('Display Segmented Image')

ylabel("Log Runtime (s)")
% 
% writematrix([timesSegment; timesSegmentNew], 'segmentTimes.csv')
% writematrix([timesExport; timesExportNew], 'exportDefaultTimes.csv')
% writematrix([timesExportAll; timesExportNewAll], 'exportAllTimes.csv')
% writematrix([timeChange; timesChangeImage], 'changeImageTimes.csv')

pSeg
pExpDef
pExpAll
pImage

%% 





