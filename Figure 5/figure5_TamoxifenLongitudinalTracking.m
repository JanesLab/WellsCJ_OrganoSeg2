% clear
% clc
% close all hidden
warning off

%% Load and organize data

cases = ["81", "85", "87", "88", "92", "93", "116a", "116b", "118a", "118b", "118c"];
casePattern = "BCO" + digitsPattern(2,3) + ("a"|"b"|"c"|"_");
wellPattern = ("B"|"C"|"D"|"E"|"F") + digitsPattern(2);
datesStandard = [2 5 7 10 12 14]; 
dates116 = [3 5 7 10 12 14];
dates = containers.Map(cases, {datesStandard, datesStandard, datesStandard, datesStandard, datesStandard, datesStandard, dates116, dates116, datesStandard, datesStandard, datesStandard });
conditionNames = ["DMSO", "-Est", "Tam Low", "Tam Med"];
conditions.case116a = containers.Map(["B02", "B03", "B04", "B05"], ["DMSO", "-Est", "Tam Low", "Tam Med"]);
conditions.case93 = containers.Map(["D03", "E02", "B02", "B04", "E05"], ["DMSO", "DMSO", "-Est", "Tam Low", "Tam Med"]);
conditions.case118b = containers.Map(["D02", "D03", "D04", "D05"], ["DMSO", "-Est", "Tam Low", "Tam Med"]);
conditions.case118a = containers.Map(["B02", "B03", "B04", "B05"], ["DMSO", "-Est", "Tam Low", "Tam Med"]);
conditions.case87 = containers.Map(["B02", "B03", "B04", "B05","B06", "B07", "B08"], ["DMSO", "DMSO", "-Est","-Est", "Tam Low","Tam Low", "Tam Med"]);
conditions.case92 = containers.Map(["B02", "C04", "B06", "B08"], ["DMSO", "-Est","Tam Low","Tam Med"]);
conditions.case88 = containers.Map(["C03", "C04", "C05","C06", "C07", "C08"], ["DMSO", "-Est","-Est", "Tam Low","Tam Low", "Tam Med"]);
conditions.case81 = containers.Map(["B02", "C02", "C03", "C04"], ["DMSO", "-Est","Tam Low","Tam Med"]);
conditions.case116b = containers.Map(["E02","E03", "E04", "E05"], ["DMSO", "-Est","Tam Low","Tam Med"]);
conditions.case118c = containers.Map(["F02", "F03", "F04", "F05"], ["DMSO", "-Est", "Tam Low", "Tam Med"]);
conditions.case85 = containers.Map(["B02", "C02", "D02", "E02"], ["DMSO", "-Est", "Tam Low", "Tam Med"]);

organoidAreaData = struct('filename',{}, 'caseNum', {}, 'wellNum', {},'organoidID',{},'I1',{}, 'I2',{},'I3',{},'I4',{}, 'I5', {},'I6',{}, 'dateSequence', {});
filenames = {dir("excelFiles").name};
nextIndex = 0;
h = waitbar(0, "Loading organoid data");
for file = 1:length(filenames)
    filename = filenames{file};
    if strcmp(filename(1), 'B')
        
        caseNum = char(erase((extract(filename, casePattern)), "BCO"));
        caseNum = erase(caseNum, "_");
        if ismember(caseNum, keys(dates))
            wellNum = char(extract(filename, wellPattern));
            dateSequence = dates(caseNum);
            condition = conditions.("case"+caseNum)(wellNum);
            data = readmatrix(fullfile("excelFiles", filename), 'Sheet', 'Area');
            data = data(2:end,:);
            for row = 1:size(data, 1)
                nextIndex = nextIndex+1;
                organoidAreaData(nextIndex).filename = filename;
                organoidAreaData(nextIndex).caseNum = caseNum;
                organoidAreaData(nextIndex).wellNum = wellNum;
                organoidAreaData(nextIndex).condition = condition;
                organoidAreaData(nextIndex).organoidID = data(row,1);
                areas = num2cell(data(row,2:end));
                [organoidAreaData(nextIndex).I1, organoidAreaData(nextIndex).I2, organoidAreaData(nextIndex).I3, organoidAreaData(nextIndex).I4, organoidAreaData(nextIndex).I5, organoidAreaData(nextIndex).I6] = areas{:}; 
                organoidAreaData(nextIndex).dateSequence = dateSequence;
            end
        end
    end
    waitbar(file/length(filenames),h)

end
close(h)

%% Subset data by condition

match_DMSO = strcmp({organoidAreaData.condition}, 'DMSO');
organoidAreaData_control = organoidAreaData(find(match_DMSO));

%% Gompertz Modeling

beta = zeros(3,length(organoidAreaData));

h = waitbar(0, "Fitting organoid traces");
for j = 1:length(organoidAreaData)
    [dates_rm,areas_rm] = getAllAreas(organoidAreaData,j);
    if length(areas_rm) > 3
        [beta(:,j),res] = lsqcurvefit(@(b,t) areagrowthmod(b,t),[(200*(areas_rm(1)/dates_rm(1))) (areas_rm(1)/dates_rm(1)) 0], ... %%QUESTION FOR KEVIN: initial values for C, A, and k??
            dates_rm, areas_rm, [0 0 0], inf(1,3), optimset('Display','off')); %Unconstrained nonlinear regression
    end
    waitbar(j/length(organoidAreaData),h)
end
close(h)




%% Visualize gompertz model fitting with heatmaps

observedData = [];
estimatedData = [];
dateSequencesFiltered = {};
indicesUsed = [];
for j = 1:length(organoidAreaData)
    t = organoidAreaData(j).dateSequence;
    C = beta(1,j);
    A = beta(2,j);
    k = beta(3,j);
    
    f=@(t) C*exp((2/3)*log(A/C)*exp(-k*t));
    
    y = feval(f,t);
    if C ~= 0 
        observedData(end+1, 1:6) = [organoidAreaData(j).I1 organoidAreaData(j).I2 organoidAreaData(j).I3 organoidAreaData(j).I4 organoidAreaData(j).I5 organoidAreaData(j).I6];
        estimatedData(end+1, 1:6) = y;
        indicesUsed(end+1) = j;
        
    end
end


%% 








 

dates_rep = [2 5 7 10 12 14];
count = 0;
for caseNum = ["118c"]%["118a", "118b", "118c", "116a", "116b", "87"]
   for conditionName = ["Tam Med"]%conditionNames
       count = count+1;
        fig = figure(1);
        match = strcmp({organoidAreaData.caseNum}, caseNum) & strcmp({organoidAreaData.condition}, conditionName);
        indicesInFiltered = [];
        relativeErrors_rep_gom = [];
        count = 0;
        for indexInOriginal = find(match)

            indexInFiltered = find(indicesUsed == indexInOriginal);
            if ~isempty(indexInFiltered)
                indexInFiltered = indexInFiltered(1);
                indicesInFiltered(end+1) = indexInFiltered;
                

                t = linspace(0,14,1000);
                C = beta(1,indexInOriginal);
                A = beta(2,indexInOriginal);
                k = beta(3,indexInOriginal);

                f=@(t) C*exp((2/3)*log(A/C)*exp(-k*t));
                y = feval(f,t)*(2.545^2);

                
                plot(t,y, 'Color', [0.5 0.5 0.5]);
                estimated = feval(f,dates_rep);
                count = count + 1;
                for day = 1:length(dates_rep)
                    if ~isnan(observedData(indexInFiltered,day))
                        relativeErrors_rep_gom(end+1) = 100*abs((observedData(indexInFiltered,day) - estimated(day))/(observedData(indexInFiltered,day)));
                    end
                end
                
                
                hold on
                title(sprintf("Case %s %s ", caseNum, conditionName))


            end
        end


        
        populationAverage = mean(observedData(indicesInFiltered, :), 'omitmissing')*(2.545^2);
        populationDev = std(observedData(indicesInFiltered, :), 'omitmissing')*(2.545^2);
        dates_rep = [2 5 7 10 12 14];
        

        [betaRepExp,res] = lsqcurvefit(@(b,t) areagrowthmod_exp(b,t),[(200*(populationAverage(1)/dates_rep(1))) (populationAverage(1)/dates_rep(1)) 0], ... 
            dates_rep, populationAverage, [0 0 0], inf(1,3), optimset('Display','off')); %Unconstrained nonlinear regression

        f_rep_Exp=@(t) betaRepExp(2)*exp(betaRepExp(3)*t);
       
        

        y_exp = feval(f_rep_Exp,t);                
        plot(t,y_exp, '--k', 'LineWidth', 3);
        

        errorbar(organoidAreaData(indexInOriginal).dateSequence, populationAverage, populationDev, '.k', 'MarkerSize', 40, 'LineWidth', 3);
        
        
        stylegraph(gca)
        pbaspect([1,1,1])
        
        hold off
        indicesToPlot = find(match);
        nonZero = find(beta(1,match)~=0);
        indicesToPlot = indicesToPlot(nonZero);
        figure(2)
        cdfplot(log10(beta(1,indicesToPlot)*(2.545^2)))
        stylegraph(gca)
        pbaspect([2,1,1])
        ylabel("CDF")
        xlabel("Carrying Capacity")
        figure(3)
        cdfplot(log10(beta(3,indicesToPlot)))
        stylegraph(gca)
        pbaspect([2,1,1])
        ylabel("CDF")
        xlabel("Growth Rate")
        

   end
end














relativeErrors = [];
for org = 1:size(observedData, 1)
    for date = 1:size(observedData, 2)
        relativeErrors(end+1) = abs((observedData(org,date) - estimatedData(org, date))/(observedData(org,date)));
    end
end
figure(4)
histogram(relativeErrors, linspace(0,3,100))
subtitle("All Data")
title("Relative error per entry")
medrelativeError = median(relativeErrors, 'omitnan');
relativeError25 = prctile(relativeErrors, 25);
relativeError75 = prctile(relativeErrors, 75);
relativeError5 = prctile(relativeErrors, 5);
relativeError95 = prctile(relativeErrors, 95);
hold on
ints{1} = xline(medrelativeError, 'black', 'LineWidth', 3);
ints{2} = xline(relativeError25, 'blue', 'LineWidth', 3);
ints{3} = xline(relativeError5, 'red', 'LineWidth', 3);
xline(relativeError75, 'blue', 'LineWidth', 3);
xline(relativeError95, 'red', 'LineWidth', 3);
legend([ints{:}], {sprintf("Med = %.2f", medrelativeError),sprintf("%.2f - %.2f", relativeError25, relativeError75),sprintf("%.2f - %.2f", relativeError5, relativeError95)})


%% 


cs = num2cell(beta(1,:)');
as = num2cell(beta(2,:)');
ks = num2cell(beta(3,:)');
[organoidAreaData.C] = cs{:};
[organoidAreaData.A] = as{:};
[organoidAreaData.k] = ks{:};







    
        
            
    
    
    
    







%% 







%% 

caseCount1 = 0;
colorsByCondition = [0 0.447 0.741;0.85 0.325 0.098; 0.920 0.694 0.125; 0.494 0.184 0.556];
 
p_s_ByCase_k = [];
p_s_ByCase_C = [];
ksTestResultsCase_ksStat.C = zeros(length(cases), length(cases));
ksTestResultsCase_ksStat.k = zeros(length(cases), length(cases));
sampleSizes = zeros(length(cases), 4);
comps = ["DMSO-Est", "-EstTam Low", "-EstTam Med", "Tam LowTam Med"];

medianGrowth = zeros(1, length(cases));
medianCC = zeros(1,length(cases));
lowQGrowth = zeros(1, length(cases));
lowQCC = zeros(1, length(cases));
highQGrowth = zeros(1, length(cases));
highQCC = zeros(1, length(cases));
numPerCase = zeros(1, length(cases));
for caseNum1 = cases
    caseCount1 = caseCount1 + 1;
    

    caseCount2 = 0;
    



    conditionName = "DMSO";
    
    
    
    match = strcmp({organoidAreaData.caseNum}, caseNum1) & strcmp({organoidAreaData.condition}, conditionName);
    subset = organoidAreaData(find(match));
    C = [subset.C];
    A = [subset.A];
    k = [subset.k];
    nonZero = find(C~=0);
    
    C = C(nonZero);
    A = A(nonZero);
    k = k(nonZero);
    numPerCase(caseCount1) = length(k);
    

    
   
 
  

    medianGrowth(caseCount1) = median(k);
    medianCC(caseCount1) = median(C);
    lowQGrowth(caseCount1) = prctile(k, 25);
    highQGrowth(caseCount1) = prctile(k, 75);
    lowQCC(caseCount1) = prctile(C, 25);
    highQCC(caseCount1) = prctile(C, 75);
    

  
    
     for caseNum2 = cases
        
        caseCount2 = caseCount2 + 1;            
        matchComp = strcmp({organoidAreaData.caseNum}, caseNum2) & strcmp({organoidAreaData.condition}, conditionName);
        subsetComp = organoidAreaData(find(matchComp));
        C_comp = [subsetComp.C];
        A_comp = [subsetComp.A];
        k_comp = [subsetComp.k];
        nonZero_comp = find(C_comp~=0);
        C_comp = C_comp(nonZero_comp);
        
        k_comp = k_comp(nonZero_comp);
        [~,p_C, ks_stat_C] = kstest2(log10(C), log10(C_comp));
        
        [~,p_k, ks_stat_k] = kstest2(log10(k), log10(k_comp));
        ksTestResultsCase_ksStat.C(caseCount1, caseCount2) = ks_stat_C;            
        ksTestResultsCase_ksStat.k(caseCount1, caseCount2) = ks_stat_k;
        if caseCount1 < caseCount2
            p_s_ByCase_k(end+1) = p_k;
            p_s_ByCase_C(end+1) = p_C;
        end
      
       



    end
   
    



end







%% 




dist = pdist(ksTestResultsCase_ksStat.k);
squareDist = squareform(dist);
Z = linkage(squareDist, 'ward');
figure(10)
[~,~, outperm] = dendrogram(Z, size(observedData,1));
xticklabels(cases(outperm))
title("Growth rate clustering")



displayk = ksTestResultsCase_ksStat.k(outperm,outperm(linspace(11,1,11)));

for case1 = 1:length(cases)
    for case2 = 1:length(cases)
        if case1 + case2 > length(cases)+1
            displayk(case1, case2) = nan;
            
        end

    end
end
figure(11)
heatmap(cases(outperm(linspace(11,1,11))),cases(outperm),round(displayk,2,"significant"));
ksTestResultsCase_ksStat.k(outperm,outperm);
colormap("gray")
clim([0,0.7])






title("Ks stat for growth rate by case")



load("populationGrowthRates.mat")
load("populationGrowthRates_low.mat")
load("populationGrowthRates_high.mat")
figure(12)
medianGrowthRearrange = medianGrowth(outperm);

lowQGrowthRearrange = lowQGrowth(outperm);

highQGrowthRearrange = highQGrowth(outperm);

dmsoBetaRearrange = betaTotal(1, outperm);
dmsoBetaRearrange_low = betalower_total(1, outperm);
dmsoBetaRearrange_high = betaupper_total(1, outperm);


x = categorical(cases(outperm));
x = reordercats(x,cases(outperm));



bar(x,dmsoBetaRearrange)
hold on
errorbar(x, dmsoBetaRearrange, dmsoBetaRearrange-dmsoBetaRearrange_low, dmsoBetaRearrange_high-dmsoBetaRearrange, '.k', 'MarkerSize', 15);
color = [0.8500 0.3250 0.0980];
hold off
figure(13)

errorbar(x, medianGrowthRearrange, medianGrowthRearrange-lowQGrowthRearrange, highQGrowthRearrange-medianGrowthRearrange, '.', 'MarkerSize', 15, 'Color', color);

stylegraph(gca)



%% 




dist = pdist(ksTestResultsCase_ksStat.C);
squareDist = squareform(dist);
Z = linkage(squareDist, 'ward');
figure(20)
[~,~, outperm] = dendrogram(Z, size(observedData,1));
xticklabels(cases(outperm))
title("Carrying capacity clustering")


displayC = ksTestResultsCase_ksStat.C(outperm,outperm(linspace(11,1,11)));
for case1 = 1:length(cases)
    for case2 = 1:length(cases)
        if case1 + case2 > length(cases)+1
            displayC(case1, case2) = nan;
            
        end

    end
end
figure(21)
heatmap(cases(outperm(linspace(11,1,11))),cases(outperm), round(displayC, 2, "significant"));
colormap("gray")
clim([0,0.7])




title("Ks stat for carrying capacity by case")
figure(23)
medianCCRearrange = medianCC(outperm)*(2.545^2);
lowQCCRearrange = lowQCC(outperm)*(2.545^2);
highQCCRearrange = highQCC(outperm)*(2.545^2);
x = categorical(cases(outperm));
x = reordercats(x,cases(outperm));

h = errorbar(x, log10(medianCCRearrange), log10(medianCCRearrange)-log10(lowQCCRearrange), log10(highQCCRearrange)-log10(medianCCRearrange), '.', 'MarkerSize', 15, 'Color', color);
stylegraph(gca)

%% 


[h_C_byCase crit_p_C_byCase adj_ci_C_byCase adj_p_C_byCase] = fdr_bh(p_s_ByCase_C, 0.05, 'pdep', 'yes');
[h_k_byCase crit_p_k_byCase adj_ci_k_byCase adj_p_k_byCase] = fdr_bh(p_s_ByCase_k, 0.05, 'pdep', 'yes');

table_adj_P_k = zeros(11);
table_adj_P_C = zeros(11);
numComparisons = 11;
startIndTable = 1;
startIndArray = 1;
for p_1 = 1:size(table_adj_P_k, 1)
    startIndTable = startIndTable+1;
    numComparisons = numComparisons - 1;
    table_adj_P_k(p_1, startIndTable:end) = adj_p_k_byCase(startIndArray:startIndArray+numComparisons-1);
    table_adj_P_C(p_1, startIndTable:end) = adj_p_C_byCase(startIndArray:startIndArray+numComparisons-1);
    startIndArray = startIndArray + numComparisons;

end






%% Define functions

% Gompertz model



function ggm = areagrowthmod(b,t)

C = b(1);
A = b(2);
k = b(3);
ggm = C*exp((2/3)*log(A/C)*exp(-k*t));

end

function ggm = areagrowthmod_exp(b,t)
A = b(2);
k = b(3);
ggm = A * exp(k*t);
end

% Area and date extraction

function [dates, areas] = getAllAreas(areaStruct, row)

    dates = areaStruct(row).dateSequence;
    areas = [areaStruct(row).I1 areaStruct(row).I2 areaStruct(row).I3 areaStruct(row).I4 areaStruct(row).I5 areaStruct(row).I6];
    combined = [dates;areas];
    combinedRemoved = rmmissing(combined, 2);
    dates = combinedRemoved(1,:);
    areas = combinedRemoved(2,:);
end

function D2 = naneucdist(XI,XJ)  
    %NANEUCDIST Euclidean distance ignoring coordinates with NaNs
    n = size(XI,2);
    sqdx = (XI-XJ).^2;
    nstar = sum(~isnan(sqdx),2); % Number of pairs that do not contain NaNs
    nstar(nstar == 0) = NaN; % To return NaN if all pairs include NaNs
    D2squared = sum(sqdx,2,'omitnan').*n./nstar; % Correction for missing coordinates
    D2 = sqrt(D2squared);
end



function stylegraph(h)
%This function changes the graph style from MATLAB defaults
    grid off %turn off grid
    h.XColor = [0 0 0]; %change x-axis color
    h.YColor = [0 0 0]; %change y-axis color
    h.TickDir = 'out'; %tick marks out
    h.TickLength = [0.02 0.05]; %increase tick length
    h.Box = 'off'; %no box
end




