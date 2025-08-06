casePattern = "BCO" + digitsPattern(3);
wellPattern = ("B"|"C"|"D") + digitsPattern(2);
reporterNames = ["Nucview", "None", "Live-Dead"];
conditionNames = ["0 Gy", "1 Gy", "6 Gy"];
imageNames = ["P1", "P2", "R1", "R2", "R3", "R4", "R5", "A1", "A2", "DAPI"];
caseNames = ["149", "157", "155", "160", "161"];
reporters.case149 = containers.Map(["C","B","D"],reporterNames);
reporters.case157 = containers.Map(["B","D","C"],reporterNames);
reporters.case155 = containers.Map(["D","B","C"],reporterNames);
reporters.case160 = containers.Map(["D","C","B"],reporterNames);
reporters.case161 = containers.Map(["B","D","C"],reporterNames);
conditions.case149 = containers.Map(["02", "03", "04"], conditionNames);
conditions.case157 = containers.Map(["04", "05", "06"], conditionNames);
conditions.case155 = containers.Map(["02", "03", "04"], conditionNames);
conditions.case160 = containers.Map(["02", "03", "04"], conditionNames);
conditions.case161 = containers.Map(["08", "09", "10"], conditionNames);
dates149_157 = [2 nan 5 6 7 8 9 12 14]; 
dates155 = [2 4 6 7 8 9 10 13 14];
dates160 = [2 4 6 7 8 9 10 12 14];
dates161 = [2 nan 5 6 7 8 9 12 14];
dates = containers.Map(caseNames, {dates149_157, dates149_157, dates155, dates160, dates161});

percentileShift = 20;

%means = mean(data, 1, 'omitnan');
% plot(1:10, means)
% boxplot(data)
organoidNucFluorData = struct('filename',{}, 'caseNum', {}, 'wellNum', {}, 'condition', {}, 'reporter',{},'organoidID',{},'P1',{}, 'P2',{},'R1',{},'R2',{}, 'R3', {},'R4',{}, 'R5', {}, 'A1', {}, 'A2', {}, 'DAPI', {}, 'dateSequence', {});
filenames = {dir("excelFiles").name};
nextIndex = 0;
h = waitbar(0, "Loading organoid data");

for file = 1:length(filenames)
    filename = filenames{file};
    if strcmp(filename(1), 'B') &contains(filename, "xls")
        
        caseNum = char(extract(extract(filename, casePattern), digitsPattern));
        wellNum = char(extract(filename, wellPattern));
        condition = conditions.("case"+caseNum)(wellNum(2:3));
        reporter = reporters.("case"+caseNum)(wellNum(1));
     
        dateSequence = dates(caseNum);
        data = readmatrix(fullfile("excelFiles", filename), 'Sheet', 'Percentile Pixel Intensity2');
        data = data(4:end,:);
        dapiData = readmatrix(fullfile("excelFiles", filename), 'Sheet', 'Percentile Pixel Intensity1');
        dapiData = dapiData(4:end,end);


        for row = 1:size(data, 1)
            nextIndex = nextIndex+1;
            organoidNucFluorData(nextIndex).filename = filename;
            organoidNucFluorData(nextIndex).caseNum = caseNum;
            organoidNucFluorData(nextIndex).wellNum = wellNum;
            organoidNucFluorData(nextIndex).condition = condition;
            organoidNucFluorData(nextIndex).reporter = reporter;
            
            organoidNucFluorData(nextIndex).organoidID = data(row,1);
            
            fluor = num2cell(data(row,2:end-2));
            if size(data,2) == 11
                
                [organoidNucFluorData(nextIndex).P1, organoidNucFluorData(nextIndex).R1, organoidNucFluorData(nextIndex).R2, organoidNucFluorData(nextIndex).R3, organoidNucFluorData(nextIndex).R4, organoidNucFluorData(nextIndex).R5, organoidNucFluorData(nextIndex).A1, organoidNucFluorData(nextIndex).A2] = fluor{:}; 
                organoidNucFluorData(nextIndex).P2 = nan;
            else             
                
               
                [organoidNucFluorData(nextIndex).P1, organoidNucFluorData(nextIndex).P2, organoidNucFluorData(nextIndex).R1, organoidNucFluorData(nextIndex).R2, organoidNucFluorData(nextIndex).R3, organoidNucFluorData(nextIndex).R4, organoidNucFluorData(nextIndex).R5, organoidNucFluorData(nextIndex).A1, organoidNucFluorData(nextIndex).A2] = fluor{:}; 
            end
            organoidNucFluorData(nextIndex).DAPI = dapiData(row);
            organoidNucFluorData(nextIndex).dateSequence = dateSequence;
            
        end

    
    end
    waitbar(file/length(filenames),h)

end
close(h)


organoidNucFluorDataPreShift = organoidNucFluorData;

for caseNum = 1:length(caseNames)
    match_Case = strcmp({organoidNucFluorDataPreShift.caseNum}, caseNames(caseNum));
    fullCaseData = organoidNucFluorDataPreShift(find(match_Case));        
    for imageNum = 1:length(imageNames)
        if ~isnan(percentileShift)
            
            fullImageDataAcrossConditions = [fullCaseData.(imageNames(imageNum))];
            shift = prctile(fullImageDataAcrossConditions, percentileShift);
            organoidsInCase = find(match_Case);
            for org = organoidsInCase
                organoidNucFluorData(org).(imageNames(imageNum)) = organoidNucFluorDataPreShift(org).(imageNames(imageNum)) - shift;
            end
            
        end
    end
       
end









%%  Gaussian Mixture fitting and histograms










addValues = .001;
binSize = 5;
gmmFigure = figure(1);
gmmTileFinal = tiledlayout(2,3);

title("Data by case gaussian mixture")
cutoffsByCase = zeros(1,length(caseNames));

for caseNum = 1:length(caseNames)
    for addValue = addValues
        match_Case = strcmp({organoidNucFluorData.caseNum}, caseNames(caseNum));
        fullCaseData = organoidNucFluorData(find(match_Case));
        
        fullCaseIntensities = [fullCaseData.R2; fullCaseData.R3; fullCaseData.R4; fullCaseData.R5; fullCaseData.A1; fullCaseData.A2]';
        
        fullCaseIntensities = fullCaseIntensities(:);
        ax = nexttile;
        binrng = min(fullCaseIntensities)-100:binSize:max(fullCaseIntensities)+100;
        counts = histc(fullCaseIntensities, binrng);

        b0 = [20 60 4 4 1 0.2]; 
        [mean1, mean2, std1, std2, mix1, mix2] = nlfitgmModel(binrng, counts, binSize, addValue, b0);
        
        x = [min(fullCaseIntensities)-100:0.1:max(fullCaseIntensities)+100];
        y1 = normpdf(x,mean1,std1)*mix1;   
        y2 = normpdf(x,mean2,std2)*mix2; 
        bar(binrng,counts/sum(counts),'histc')
        hold on
        plot(x, y1, 'LineWidth', 2) 

        plot(x, y2, 'LineWidth', 2)
        xline(mean1+std1*1.65, 'LineWidth', 2, 'Color', [0,0,0])
        cutoffsByCase(caseNum) = mean1+std1*1.645;

        ylim([0 max(counts/sum(counts))])
        xlim([-100 500])
        legend({'', sprintf('\\mu1 = %d, \\sigma1 = %d', round(mean1), round(std1)), sprintf('\\mu2 = %d, \\sigma2 = %d', round(mean2), round(std2)), sprintf('Cutoff = %.1f',mean1+std1*1.645)})

        pbaspect(ax, [1,1,1])
        title(sprintf("Case %s, All Conditions Together", caseNames(caseNum)))


        hold off

        
        
    end
end
hold off






%% Get DAPI Cutoffs

cutoffsByCaseNucDapi = zeros(1,length(caseNames));
figure(2)
tiledlayout(2,3)
addValue = .005;
binSize = 10;



    
for caseNum = 1:length(caseNames)
   
    match_Case = strcmp({organoidNucFluorData.caseNum}, caseNames(caseNum));
    fullCaseData = organoidNucFluorData(find(match_Case));
    
    dapiData = [fullCaseData.DAPI]';
   
    
    binrng = min(dapiData)-500:binSize:max(dapiData)+500;  
    counts = histc(dapiData, binrng);
    
        [mean1, mean2, std1, std2, mix1, mix2] = nlfitgmModelDAPI(binrng, counts, binSize, addValue);
        
        x = [min(dapiData)-500:.1:max(dapiData)+500];
        y1 = normpdf(x,mean1,std1)*mix1;   
        y2 = normpdf(x,mean2,std2)*mix2; 
        nexttile
        bar(binrng, counts/sum(counts),'histc')
        
        hold on

        plot(x, y1, 'LineWidth', 2) 

        plot(x, y2, 'LineWidth', 2)
        
        xline(mean1+std1*1.645, 'LineWidth', 2)
       
        cutoffsByCaseNucDapi(caseNum) = mean1+std1*1.645;
        
        
        xlim([-400 1500])

        title(sprintf("Case Number: %s", caseNames(caseNum)))
        legend({'', sprintf('\\mu1 = %d, \\sigma1 = %d', round(mean1), round(std1)), sprintf('\\mu2 = %d, \\sigma2 = %d', round(mean2), round(std2)), sprintf('Cutoff = %.1f',mean1+std1*1.645)})
    
    
    
    
    hold off
    
    
    
end






%% Individual Organoid Traces


numOrganoidsTracked = zeros(length(caseNames), length(conditionNames));
numOrganoidsTrackedTrimmed = zeros(length(caseNames), length(conditionNames));





figure(3)
percentDeadByCaseAndCondition = zeros(length(caseNames), length(conditionNames));
percentDeadByTimeAndCondition = cell(1,length(caseNames));
predictorAndTime = cell(1,length(caseNames));
organoidTraces = tiledlayout(length(caseNames), length(conditionNames));
flippedOrganoids = {};
map = [0 0 0; .75, .75, .75; 0 1 0; 0 0 1];
fullOrganoidMatrix = [];
fullFirstDeadImage = [];




count = 0;
doImageCount = 0;
for caseNum = 1:5
    percentDeadTimeConditionThisCase = zeros(length(conditionNames), 6);
    predictorAndTimeThisCase = [];
    for condition = 1:length(conditionNames)
        match_CaseAndCondition = strcmp({organoidNucFluorData.caseNum}, caseNames(caseNum)) & strcmp({organoidNucFluorData.condition}, conditionNames(condition));
        fullCaseAndConditionData = organoidNucFluorData(find(match_CaseAndCondition));
        fullCaseAndConditionIntensities = [fullCaseAndConditionData.R2; fullCaseAndConditionData.R3; fullCaseAndConditionData.R4; fullCaseAndConditionData.R5; fullCaseAndConditionData.A1; fullCaseAndConditionData.A2]';
        groupNumOrg = size(fullCaseAndConditionIntensities, 1);
        organoidMatrix = zeros(groupNumOrg, size(fullCaseAndConditionIntensities, 2)+1);
        firstDeadImageByOrg = zeros(1,length(groupNumOrg));
        fullCaseAndConditionIntensitiesDAPI = [fullCaseAndConditionData.DAPI]';

        for organoid = 1:groupNumOrg


            firstDeadImage = nan;
            for image = 1:size(fullCaseAndConditionIntensities, 2)
                if isnan(fullCaseAndConditionIntensities(organoid, image))

                    organoidMatrix(organoid, image) = 0.5;

                end
                if fullCaseAndConditionIntensities(organoid, image) >= cutoffsByCase(caseNum)
                    organoidMatrix(organoid, image:end-1) = 1;
                    if isnan(firstDeadImage)
                        firstDeadImage = image;

                    end
                end

            end

            if isnan(fullCaseAndConditionIntensitiesDAPI(organoid))
                organoidMatrix(organoid, end) = 0.5;

            elseif fullCaseAndConditionIntensitiesDAPI(organoid) >= cutoffsByCaseNucDapi(caseNum)
                organoidMatrix(organoid, end) = 1.5;


                if isnan(firstDeadImage)
                    firstDeadImage = image+1;

                end
            end

            if isnan(firstDeadImage)
                firstDeadImage = image+2;

            end
            firstDeadImageByOrg(organoid) = firstDeadImage;
            fullFirstDeadImage(end+1) = firstDeadImage;

        end

        [sorted, index] = sort(firstDeadImageByOrg);
        organoidMatrixSorted = organoidMatrix(index,:);
        numNan = zeros(1,size(organoidMatrixSorted,1));
        for row = 1:size(organoidMatrixSorted,1)
            numNan(row) = length(find(organoidMatrixSorted(row,:) == 0.5));
        end

        organoidMatrixSortedTrimmed = organoidMatrixSorted(numNan<4,:);
        numDead = sum(organoidMatrixSortedTrimmed == 1);
        numAlive = sum(organoidMatrixSortedTrimmed == 0);
        percentDeadByCaseAndCondition(caseNum, condition) = numDead/(numAlive + numDead) * 100;
        for timePoint = 1:6
            numDeadTime = sum(organoidMatrixSortedTrimmed(:,timePoint) == 1);
            numAliveTime = sum(organoidMatrixSortedTrimmed(:,timePoint) == 0);
            percentDeadTimeConditionThisCase(condition, timePoint) = numDeadTime/(numAliveTime + numDeadTime) * 100;

        end

        for org = 1:size(organoidMatrixSortedTrimmed, 1)
            firstDead = find(organoidMatrixSortedTrimmed(org, :) == 1, 1,'first');
            if ~isempty(firstDead)
                date = fullCaseAndConditionData(1).dateSequence(firstDead+3);
                predictorAndTimeThisCase(end+1, :) = [condition, date, 0];

            else
                predictorAndTimeThisCase(end+1, :) = [condition, 14, 1]; %for organoids that don't die (nan or censored value)
            end


        end
        makeImages = false;
        imageNum = 0;
        if condition == 1 & caseNum == 1
            try 
                load("BCO149_C02_Fluorescence_4.mat")
                makeImages = true;
                imageNum = 11;
            catch
                "Mat files not in folder, not displaying representative images"
                
            end
               

        elseif condition == 2 & caseNum == 1
            try
                load("BCO149_C03_Fluorescence_4.mat")
                makeImages = true;
                imageNum = 12;
            catch
                "Mat files not in folder, not displaying representative images"
                
            end
            


        end
        if makeImages

            lastNuc = hiddenProps.f_displayBulkI{end-7}; 
        
            lastNuc = double(lastNuc);
            scaleHigh = prctile(lastNuc(:), 99.5);
            scaleLow = prctile(lastNuc(:), 1);

            


            dapi = hiddenProps.f_displayBulkI{end-2};
            dapi = double(dapi);
            dapi = (dapi - prctile(dapi(:),1))./(prctile(dapi(:), 99) - prctile(dapi(:),1));

            ccNuc = hiddenProps.bf_cc{end - 2};
            ccDapi = hiddenProps.bf_cc{end};
            goodOrgNumbers = [fullCaseAndConditionData.organoidID];
            nullNuc = goodOrgNumbers(find(organoidMatrix(:,end-1)==0.5));
            negativeNuc = goodOrgNumbers(find(organoidMatrix(:,end-1)==0));
            nullDapi = goodOrgNumbers(find(organoidMatrix(:,end)==0.5));
            negativeDapi = goodOrgNumbers(find(organoidMatrix(:,end)==0));
            bothPos = goodOrgNumbers(find(organoidMatrix(:,end)==1.5 & organoidMatrix(:,end-1)==1));
            matchingIDsNuc = hiddenProps.bf_matchingID{end-2};
            matchingIDsDapi = hiddenProps.bf_matchingID{end};


            ccPosNuc = ccNuc;
            ccNegNuc = ccNuc;
            ccBothNuc = ccNuc;
            ccDapiNuc = ccNuc;
            ccPosDapi = ccDapi;
            ccNegDapi = ccDapi;
            ccBothDapi = ccDapi;
            for orgNuc = 1:length(ccNuc.PixelIdxList)
                if ismember(matchingIDsNuc(orgNuc), nullNuc)
                    ccPosNuc.PixelIdxList{orgNuc} = [];
                    ccNegNuc.PixelIdxList{orgNuc} = [];
                end
                if ismember(matchingIDsNuc(orgNuc), negativeNuc)
                    ccPosNuc.PixelIdxList{orgNuc} = [];
                end
                if ~ismember(matchingIDsNuc(orgNuc), bothPos)
                    ccBothNuc.PixelIdxList{orgNuc} = [];
                end
                if ismember(matchingIDsNuc(orgNuc), nullDapi)
                    ccDapiNuc.PixelIdxList{orgNuc} = [];

                end
                if ismember(matchingIDsNuc(orgNuc), negativeDapi)
                    ccDapiNuc.PixelIdxList{orgNuc} = [];
                end

            end

            for orgDapi = 1:length(ccDapi.PixelIdxList)
                if ismember(matchingIDsDapi(orgDapi), nullDapi)
                    ccPosDapi.PixelIdxList{orgDapi} = [];
                    ccNegDapi.PixelIdxList{orgDapi} = [];
                end
                if ismember(matchingIDsDapi(orgDapi), negativeDapi)
                    ccPosDapi.PixelIdxList{orgDapi} = [];
                end
                if ~ismember(matchingIDsDapi(orgDapi), bothPos)
                    ccBothDapi.PixelIdxList{orgDapi} = [];
                end
            end

            

            
            

            figure(imageNum)
            repImage = tiledlayout(2,4);
            if condition == 1
                title(repImage, "Figure6b - Case 149 0Gy")
            else
                title(repImage, "Figure6c - Case 149 1Gy")
            end



            for nucNum = 1:6
                
                nucNumProps = nucNum + 2;
                date = organoidNucFluorData(1).dateSequence(nucNumProps+1);
                nucInd = 2+(nucNumProps-1)*3;
                nucImage = hiddenProps.f_displayBulkI{nucInd};
                nucImage = double(nucImage);
                nucImage = (nucImage - scaleLow)./(scaleHigh - scaleLow);
                ccNuc = hiddenProps.bf_cc{nucNumProps};

                goodOrgNumbers = [fullCaseAndConditionData.organoidID];
                nullNuc = goodOrgNumbers(find(organoidMatrix(:,nucNum)==0.5));
                negativeNuc = goodOrgNumbers(find(organoidMatrix(:,nucNum)==0));
                bothPos = goodOrgNumbers(find(organoidMatrix(:,end)==1.5 & organoidMatrix(:,nucNum)==1));

                matchingIDsNuc = hiddenProps.bf_matchingID{nucNumProps};



                ccPosNuc = ccNuc;
                ccNegNuc = ccNuc;
                ccBothNuc = ccNuc;



                for orgNuc = 1:length(ccNuc.PixelIdxList)
                    if ismember(matchingIDsNuc(orgNuc), nullNuc)
                        ccPosNuc.PixelIdxList{orgNuc} = [];
                        ccNegNuc.PixelIdxList{orgNuc} = [];
                    end
                    if ismember(matchingIDsNuc(orgNuc), negativeNuc)
                        ccPosNuc.PixelIdxList{orgNuc} = [];
                    end
                    if ~ismember(matchingIDsNuc(orgNuc), bothPos)
                        ccBothNuc.PixelIdxList{orgNuc} = [];
                    end



                end


                nexttile
                
                allOverlayNuc = imoverlay(nucImage, imdilate(bwperim(labelmatrix(ccNegNuc)), ones(3)), [1 1 1]);
                positiveOverlayNuc = imoverlay(allOverlayNuc, imdilate(bwperim(labelmatrix(ccPosNuc)), ones(6)), [0 1 0]);
                finalOverlay = imoverlay(positiveOverlayNuc, imdilate(bwperim(labelmatrix(ccBothNuc)), ones(6)), [0 1 1]);
                imshow(finalOverlay)
                title("Day" + num2str(date))







            end

            nexttile
            
            allOverlayDapi = imoverlay(dapi, imdilate(bwperim(labelmatrix(ccNegDapi)), ones(3)), [1 1 1]);
            positiveOverlayDapi = imoverlay(allOverlayDapi, imdilate(bwperim(labelmatrix(ccPosDapi)), ones(6)), [0 0 1]);
            finalOverlay = imoverlay(positiveOverlayDapi, imdilate(bwperim(labelmatrix(ccBothDapi)), ones(6)), [0 1 1]);
            imshow(finalOverlay)
            title("DAPI")

        end
        figure(3)



        fullOrganoidMatrix = [fullOrganoidMatrix; organoidMatrix];



        numOrganoidsTrackedTrimmed(caseNum, condition) = size(organoidMatrixSortedTrimmed,1);

        ax = nexttile;
        imagesc(organoidMatrixSortedTrimmed)
        colormap(map)
        xticks([1 2 3 4 5 6 7])
        %pbaspect(ax, [1,1,1])
        dates = fullCaseAndConditionData(1).dateSequence;
        labelDates = [dates(4:end) "DAPI"];
        xticklabels(labelDates)
        xlabel("Days in culture")
        yticks('')
        ylabel('Individual Organoids')
        %pbaspect(ax, [1,1,1])
        title(sprintf("Case # %s, %s", caseNames(caseNum), conditionNames{condition}))

    end
    percentDeadByTimeAndCondition{caseNum} = percentDeadTimeConditionThisCase;
    predictorAndTime{caseNum} = predictorAndTimeThisCase;
end







[sorted, index] = sort(fullFirstDeadImage);
%%%%Original data no trimming, all caseNames/conditions, flipped organoids
%%%% filled in
organoidMatrixSorted = fullOrganoidMatrix(index,:);


lastRowBeforeDapi = 0;
dapiOverlap = 0;
for row = 1:size(organoidMatrixSorted,1)
    
    if organoidMatrixSorted(row, 6) == 1
        lastRowBeforeDapi = row;
    end
    
end


for row = 1:lastRowBeforeDapi
    
    if organoidMatrixSorted(row, 7) == 1.5
        dapiOverlap = dapiOverlap + 1;
    end
    
end


%% 



%% 

numNan = zeros(1,size(organoidMatrixSorted,1));
for row = 1:size(organoidMatrixSorted,1)
    numNan(row) = length(find(organoidMatrixSorted(row,:) == 0.5));
end

organoidMatrixSortedTrimmed = organoidMatrixSorted(numNan<4,:);
%%%% >50% NaN trimmed, all cases/conditions, flipped organoids
%%%% filled in
figure(4)
imagesc(organoidMatrixSortedTrimmed)
colormap(map)
xticks([1 2 3 4 5 6 7])
title("Organoid Traces All Cases With DAPI")
subtitle("Remove rows with >50% Nan")
lastRowBeforeDapi = 0;
dapiOverlap = 0;
for row = 1:size(organoidMatrixSortedTrimmed,1)
    
    if organoidMatrixSortedTrimmed(row, 6) == 1
        lastRowBeforeDapi = row;
    end
    
end


for row = 1:lastRowBeforeDapi
    
    if organoidMatrixSortedTrimmed(row, 7) == 1.5
        dapiOverlap = dapiOverlap + 1;
    end
    
end
%% 



numDapi = length(find(organoidMatrixSortedTrimmed(:,7)==1.5));
p = hygecdf(dapiOverlap, size(organoidMatrixSortedTrimmed,1), lastRowBeforeDapi, numDapi, 'upper')


%% 



cxmdl = cell(1,length(caseNames));

for i = 1:length(caseNames)
    test = cell2mat(predictorAndTime(1,i));
    test = array2table(test);
    test = renamevars(test,"test1","Predictors");
    test = renamevars(test,"test2","Times");
    test = renamevars(test,"test3","CensoredStatus");
    test.Predictors = categorical(test.Predictors);
    cxmdl{1,i} = fitcox(test,"Times",'Censoring',"CensoredStatus");
end

figure(40)
for i = 1:length (caseNames)
    %figure(40+i)
    subplot(2,3,i);
    test = cell2mat(predictorAndTime(1,i));
    test = array2table(test);
    test = renamevars(test,"test1","Predictors");
    test = renamevars(test,"test2","Times");
    test = renamevars(test,"test3","CensoredStatus");
    test0Gy = test.Predictors == 1;
    j = sum(test0Gy);
    test1Gy = test.Predictors == 2;
    k = sum(test1Gy);
    test6Gy = test.Predictors == 3;
    m = sum(test6Gy);
    ecdf(test.Times(1:j),'censoring',test.CensoredStatus(1:j),'function','survivor');
    hold on
    ecdf(test.Times(j+1:j+k),'censoring',test.CensoredStatus(j+1:j+k),'function','survivor');
    hold on
    ecdf(test.Times(j+k+1:j+k+m),'censoring',test.CensoredStatus(j+k+1:j+k+m),'function','survivor');
    
    leg = {'0Gy' sprintf('1Gy -  %.3f', cxmdl{1,i}.Coefficients{"Predictors_2","pValue"} ) sprintf('6Gy -  %.3f', cxmdl{1,i}.Coefficients{"Predictors_3","pValue"} )};
    legend(leg,'Location','SouthWest');
    title(caseNames(i));
    pbaspect([1 1 1]);
    
    
    stylegraph(gca)
    ylabel("Survival Fraction")
    xlabel("Days in Culture")
    hold off
    ylim([0,1])
    xlim([5,14])
    
end



function [mean1, mean2, std1, std2, mix1, mix2] = nlfitgmModel(bins, counts, binSize, addValue, b0)
    normModel = @(b,x) b(5)*normpdf(x,b(1),b(3))+b(6)*normpdf(x,b(2),b(4))+addValue;
    
    
    
    countsFrac = (counts/(sum(counts)))'+addValue;
   
   
    beta = nlinfit(bins, countsFrac, normModel, b0, 'ErrorModel', 'proportional');
    mean1 = beta(1);
    mean2 = abs(beta(2));
    std1 = beta(3);
    std2 = beta(4);
    mix1 = beta(5);
    mix2 = beta(6);

end

function [mean1, mean2, std1, std2, mix1, mix2] = nlfitgmModelDAPI(bins, counts, binSize, addValue)
    normModel = @(b,x) b(5)*normpdf(x,b(1),b(3))+b(6)*normpdf(x,b(2),b(4))+addValue;
    b0 = [10 200 50 50 10 0]; 
  
    countsFrac = (counts/(sum(counts)))'+addValue;
    
    sum1 = sum(countsFrac);
    
    beta = nlinfit(bins, countsFrac, normModel, b0, 'ErrorModel', 'proportional');
    mean1 = beta(1);
    mean2 = abs(beta(2));
    std1 = beta(3);
    std2 = beta(4);
    mix1 = beta(5);
    mix2 = beta(6);

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
