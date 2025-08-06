
parallel = [3.66 4.36 3.88 6.08 3.86 4.65 4.73 4.50 3.73 3.47;
    6.58 12.79 8.39 7.38 6.44 5.54 5.12 4.14 6.04 7.06;
    10.21 14.72 14.42 15.15 11.51 13.06 14.93 12.70 13.81 13.61;
    22.27 22.16 24.40 24.06 25.27 22.20 22.59 24.31 23.14 23.54;
    50.54 46.66 53.74 46.82 48.68 45.27 48.55 46.56 50.34 51.15;
    97.16 97.70 101.55 102.89 109.27 97.65 114.91 107.21 111.31 111.38];


org2 = [3.50 3.63 3.94 3.65 4.11 3.85 4.15 3.60 3.64 4.04;
    4.44 10.15 8.41 8.90 7.89 8.47 7.77 6.66 7.60 7.01;
    15.72 19.26 18.45 19.70 17.77 18.25 18.55 16.83 19.38 18.85;
    37.51 37.69 37.16 36.89 35.93 35.88 32.14 36.26 36.30 35.20;
    98.42 89.66 112.50 94.95 91.74 92.40 102.78 91.95 95.16 97.83;
    197.83 212.21 200.04 203.02 199.45 189.05 224.72 212.75 215.20 207.08];

org1 = [1.91 2.02 1.61 1.86 1.76 1.53 1.51 1.98 1.80 1.80;
    2.32 5.06 4.25 3.76 3.85 3.24 3.11 3.50 2.94 3.68;
    8.68 13.92 13.68 15.37 9.70 8.53 13.75 9.10 9.70 8.05;
    18.21 25.57 28.01 25.01 23.01 18.61 24.48 38.82 28.22 27.74;
    54.64 63.33 76.55 68.90 47.41 63.10 58.81 63.14 56.26 57.37;
    119.88 122.98 108.83 122.59 134.06 113.37 121.21 112.90 119.25 112.85];





figure(10)
mean(org1, 2)
std(org1,0,2)
h =  errorbar([1,2,5, 10, 25, 50], mean(org1,2), std(org1,0,2), '-o');
h.MarkerFaceColor = h.Color;
hold on
h = errorbar([1,2,5, 10, 25, 50], mean(org2,2), std(org2,0,2), '-o');
h.MarkerFaceColor = h.Color;
h = errorbar([1,2,5, 10, 25, 50], mean(parallel,2), std(parallel,0,2), '-o');
h.MarkerFaceColor = h.Color;
xlim([0,55])
hold off
legend("OrganoSeg", "OrganoSeg2", "OrganoSeg2 With Parallel Threading")
stylegraph(gca)
xlabel("Number of Images")
ylabel("Segmentation Runtime (s)")
pbaspect([1,1,1])
hold off

dataSetSizes = [1 2 5 10 25 50];
dataSetSizeCol = zeros(10*length(dataSetSizes),1);
organoSegCol = zeros(10*length(dataSetSizes),1);
organoSeg2ParCol = zeros(10*length(dataSetSizes),1);
organoSeg2Col = zeros(10*length(dataSetSizes),1);

for n = 1:length(dataSetSizes)
    dataSetSizeCol(1+10*(n-1): 10*n) = repmat(dataSetSizes(n), 10, 1);
    organoSegCol(1+10*(n-1): 10*n) = org1(n,:)';
    organoSeg2Col(1+10*(n-1): 10*n) = org2(n,:)';
    organoSeg2ParCol(1+10*(n-1): 10*n) = parallel(n,:)';


end
%% 

t = table(dataSetSizeCol,organoSegCol, organoSeg2ParCol, organoSeg2Col,...
'VariableNames',{'DataSetSize','OrganoSeg','OrganoSeg2Parallel','OrganoSeg2'});

platforms = table(categorical(["OrganoSeg", "OrganoSeg2Parallel", "OrganoSeg2"])',VariableNames="Platform");
rm = fitrm(t,"OrganoSeg-OrganoSeg2~DataSetSize",'WithinDesign',platforms)
ranovatbl = ranova(rm)
multcompare(rm, 'Platform')

function stylegraph(h)
%This function changes the graph style from MATLAB defaults
    grid off %turn off grid
    h.XColor = [0 0 0]; %change x-axis color
    h.YColor = [0 0 0]; %change y-axis color
    h.TickDir = 'out'; %tick marks out
    h.TickLength = [0.02 0.05]; %increase tick length
    h.Box = 'off'; %no box
end