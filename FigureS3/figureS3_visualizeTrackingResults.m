cases = [116, 116, 118, 118, 81, 85, 87, 87, 92, 93];
wellNums = ["B02", "E02", "B02", "D02", "C04", "B02", "B02", "B03","B02","D03"];
registrations = ["Rigid", "Non-rigid"];
iterations = ["Rigid", "RigidRotate", "NonRigid_Smooth05", "NonRigid_Smooth1", "NonRigid_Smooth2", "NonRigid_Smooth3"];
iterationsNames = ["Rigid", "Rigid Rotate", "NonRigid Smooth = 0.5", "NonRigid Smooth = 1", "NonRigid Smooth = 2", "NonRigid Smooth = 3"];
imageMethods = ["Centroids", "Image", "Mask"];
trackingAnalysisStructArray = struct("imageSet", {}, "imageMethod", {}, "registration", {}, "percentCorrectComparisons", {}, "percentRetention", {});
colors = [1 0 0; 0 0 0; 0 1 0; 0 0 1; 1 0 1; 1 0.6 0];
markers = ["o", "*", " +"];
clf reset
load("trackingAnalysisResults.mat")
%Uncomment for case-by-case results
for imageSet = 1:length(wellNums)
    % figure(imageSet)
    % tiledlayout(2,3)
    for imageMethodNum = 1:length(imageMethods)
        imageMethod = imageMethods(imageMethodNum);
        for iterationNum = 1:length(iterations)
            iteration = iterations(iterationNum);
            thisIterationResults = trackingResults{imageSet}.(imageMethod).(iteration);
            trackingAnalysisStructArray(end+1).imageSet = num2str(cases(imageSet)) + wellNums(imageSet);
            trackingAnalysisStructArray(end).imageMethod = imageMethod;
            trackingAnalysisStructArray(end).registration = iteration;
            percentComparisons = 100*thisIterationResults.correctComparisons./thisIterationResults.totalComparisons;
            trackingAnalysisStructArray(end).percentComparisons = percentComparisons(2:6);
            trackingAnalysisStructArray(end).percentRetained = 100*thisIterationResults.correctRetained./thisIterationResults.totalRetained;
            % nexttile(imageMethodNum)
            % plot(2:6, trackingAnalysisStructArray(end).percentComparisons, 'color', colors(iterationNum, :), 'Marker', 'o')
            
            % hold on
            % title("Percent Correct Comparisons")
            % subtitle(imageMethod)
            % ylim([0 100])
            % pbaspect([1,1,1])
            % stylegraph(gca)
            % 
            % nexttile(imageMethodNum + 3)
            % plot(1:6, trackingAnalysisStructArray(end).percentRetained, 'color', colors(iterationNum, :), 'Marker', 'o')
            % hold on
            % title("Percent Retention")
            % 
            % subtitle(imageMethod)
            % ylim([0 100])
            % pbaspect([1,1,1])
            % stylegraph(gca)

            matchingID = trackingResults{imageSet}.(imageMethod).(iteration).matchingIDs;
            maxID = 0;
            for day = 1:length(matchingID)
                if max(matchingID{day}) > maxID
                    maxID = max(matchingID{day});
                end
            end
           










        end
    end
    % for i = 1:length(colors)
    %     legendPlots{i} = plot(nan, 'color', colors(i,:));
    %     legendNames{i} = iterationsNames(i);
    % end
    % 
    % 
    % legend([legendPlots{:}], legendNames)
end






%% 
figure(12)
iterationNamesHere = ["Rigid No Rotation", "Rigid With Rotation", "Non-Rigid"];
ax = tiledlayout(2,3);
selectedIterations = [1,2,6];
methodsColors = [27,158,119;
217,95,2;
117,112,179]/255;
for iterationNum = 1:length(selectedIterations)
    selectedIterationNum = selectedIterations(iterationNum);
    iteration = iterations(selectedIterationNum);
    for imageMethodNum = 1:length(imageMethods)
        imageMethod = imageMethods(imageMethodNum);


        match = strcmp([trackingAnalysisStructArray.registration], iteration) & strcmp([trackingAnalysisStructArray.imageMethod], imageMethod);
        subset = trackingAnalysisStructArray(match);

        percentComparisonsTotal = cell2mat({subset.percentComparisons}');
        percentRetentionTotal = cell2mat({subset.percentRetained}');

        ax = nexttile(iterationNum);
        errorbar(2:6, mean(percentComparisonsTotal), std(percentComparisonsTotal)/sqrt(size(percentComparisonsTotal, 1)), 'color', methodsColors(imageMethodNum, :), 'Marker', 'o')
        hold on
        title("Average Correct Comparisons")
        subtitle(iterationNamesHere(iterationNum))
        ylim([0 100])
        pbaspect([1,1,1])


        nexttile(iterationNum + 3)
        errorbar(1:6, mean(percentRetentionTotal), std(percentRetentionTotal)/sqrt(size(percentComparisonsTotal, 1)), 'color', methodsColors(imageMethodNum, :), 'Marker', 'o')
        hold on
        title("Average Retention")
        subtitle(iterationNamesHere(iterationNum))
        ylim([0 100])
        pbaspect([1,1,1])


    end


end
legendPlots = {};
legendNames = {};
for i = 1:length(methodsColors)
    legendPlots{i} = plot(nan, 'color', methodsColors(i,:));
    legendNames{i} = imageMethods(i);
end
legend([legendPlots{:}], legendNames)


%% 

crops = {[386,1286,571,1471], [1,901, 553, 1453], [318,1218, 780,1680]};
images = [1,4,6];
for i = 1:length(images)
   
    visualizeTracings(trackingResults, images(i), wellNums, cases, "Mask", "Rigid", crops{i})
end
%% 







function visualizeTracings(trackingResults, wellInd, wellNums,caseNums, imageMethod, iteration, crop)
wellNum = wellNums(wellInd);
caseNum = num2str(caseNums(wellInd));
filenames = {dir().name};


fileInd = contains(filenames, wellNum) & contains(filenames, "BCO"+caseNum) & contains(filenames, "Tracking");
try
    loadedProps = load(filenames{fileInd}, 'hiddenProps');
    hiddenProps = loadedProps.hiddenProps;
    
    matchingIDs = trackingResults{wellInd}.(imageMethod).(iteration).matchingIDs;
    combinedTraces = trackingResults{wellInd}.combinedTraces;
    automatedTraceWithGlobalID = zeros(size(combinedTraces));
    automatedTraceWithLocalID = zeros(size(combinedTraces));
    firstAppearanceDay = zeros(1, size(combinedTraces, 1));
    organoidFirstDayIDs = zeros(1, size(combinedTraces, 1));
    
    
    
    
    for i = 1:size(combinedTraces, 2)
        
        for j = 1:length(matchingIDs{i})
            
            if ismember(j,combinedTraces(:,i))            
                traceNumber = find(combinedTraces(:,i) == j, 1, 'first');
                traceNumber = traceNumber(1);
                automatedTraceWithGlobalID(traceNumber, i) = matchingIDs{i}(j);            
            end
    
           
            
        end
    
    end
    
    
    
    for org = 1:size(combinedTraces, 1)
        firstAppearance= find(~isnan(combinedTraces(org,:)), 1, 'first');
        firstAppearanceDay(org) = firstAppearance(1);
        organoidFirstDayIDs(org) = automatedTraceWithGlobalID(org, firstAppearance(1));
    end
    
    for i = 1:size(combinedTraces, 2)
        
        for j = 1:length(matchingIDs{i})
    
            if ismember(matchingIDs{i}(j), organoidFirstDayIDs)
                traceNumber = find(organoidFirstDayIDs == matchingIDs{i}(j), 1, 'first');
                traceNumber = traceNumber(1);
                if firstAppearanceDay(traceNumber) <= i
                    automatedTraceWithLocalID(traceNumber, i) = j;    
                end
    
            end
            
        end
    
    end
    
    
    figure(20+wellInd)
    tiledlayout(2,3)
    for date = 1:size(combinedTraces, 2)
        image = hiddenProps.bulkI{date};
        image = double(image);
        image = (image - min(image(:)))./(max(image(:)) - min(image(:)));
        cc = trackingResults{wellInd}.cc{date};
        statsBW = regionprops(cc, 'Centroid');
        centersOrig = [statsBW.Centroid];
        centers = reshape(centersOrig, 2, length(statsBW))';
        basePerim = getPerimFromCC(cc);
        overlay = imoverlay(image, basePerim, [1,1,1]);
        goodOrgs = [];
        notMatchedOrgs = [];
        misMatchedOrgs = [];
        for org = 1:size(combinedTraces, 1)
            firstDay = firstAppearanceDay(org);
            if firstDay == date
                goodOrgs(end+1) = combinedTraces(org,date);
            elseif isnan(combinedTraces(org, date))
                if automatedTraceWithLocalID(org,date)~=0
                    misMatchedOrgs(end+1) = automatedTraceWithLocalID(org,date);
                end
            else
                if automatedTraceWithGlobalID(org,date) == automatedTraceWithGlobalID(org, firstDay)
                    goodOrgs(end+1) = combinedTraces(org,date);
    
                else
                    notMatchedOrgs(end+1) = combinedTraces(org,date);
                    if automatedTraceWithLocalID(org,date) ~= 0
                        misMatchedOrgs(end+1) = automatedTraceWithLocalID(org,date);
                    end
    
                end
            end
    
    
    
        end
        if ~isempty(goodOrgs)
    
            goodOrgsCC = cc;
            goodOrgsCC.PixelIdxList = cc.PixelIdxList(goodOrgs);
            goodOrgsCC.NumObjects = length(goodOrgsCC.PixelIdxList);
    
            overlay = imoverlay(overlay, getPerimFromCC(goodOrgsCC), [0,0,1]);
        end
        if ~isempty(notMatchedOrgs) 
    
            notMatchedOrgsCC = cc;
            notMatchedOrgsCC.PixelIdxList = cc.PixelIdxList(notMatchedOrgs);
            notMatchedOrgsCC.NumObjects = length(notMatchedOrgsCC.PixelIdxList);
            overlay = imoverlay(overlay, getPerimFromCC(notMatchedOrgsCC), [1,0,1]);
        end
        if ~isempty(misMatchedOrgs)
        
            misMatchedOrgsCC = cc;
            misMatchedOrgsCC.PixelIdxList = cc.PixelIdxList(misMatchedOrgs);
            misMatchedOrgsCC.NumObjects = length(misMatchedOrgsCC.PixelIdxList);
            overlay = imoverlay(overlay, getPerimFromCC(misMatchedOrgsCC), [1,1,0]);
        end
        nexttile
       
        imshow(overlay(crop(1):crop(2),crop(3):crop(4), :))
      
        
    
    
    
    
    end
catch
    "Mat files not in folder, not displaying representative images"
end




end

function perim = getPerimFromCC(cc)
label = labelmatrix(cc);
perim = bwperim(logical(label));
perim = imdilate(perim, ones(3));
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
