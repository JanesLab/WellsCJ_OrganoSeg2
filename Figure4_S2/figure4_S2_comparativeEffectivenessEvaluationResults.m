imageData = struct('filename',{}, 'imageName', {}, 'imageSource', {}, 'imageSegmenter', {}, 'trained', {}, 'image',{}, 'IOU_CB', {}, 'IOU_indivOrganoids_CB', {}, 'TP', {}, 'FP', {}, 'minArea', {}, 'Area', {});
mainDir = dir;
imageCount = 0;
colors7 = [127 201 127; 48 213 200; 255 255 40; 191 91 23;  240,2,127; 191 135 212; 255, 194, 14]/255;
%% 
h = waitbar(0, "Loading");
for i = 1:length(mainDir)
    
    if contains(mainDir(i).name, 'Images')         
        cd(mainDir(i).name)        
        cd images
        imageDir = dir;
        imageDir(1).folder
        h = waitbar(0, h,"Loading Raw Images");
        for image = 1:length(imageDir)
            
           
            if ~imageDir(image).isdir
                filename = imageDir(image).name;
                if ~contains(filename, "MATLAB")

                
                    imageCount = imageCount+1;
                    
                    imageData(imageCount).filename = filename;
                    imageName = split(imageDir(image).name, '.');
                    imageData(imageCount).imageName = imageName{1};
                    if contains(imageDir(image).folder, '\')
                        folder = split(imageDir(image).folder, '\');
                    elseif contains(imageDir(image).folder, '/')
                        folder = split(imageDir(image).folder, '/');
                    end
                    imageData(imageCount).imageSource = folder{end-1};
                    imageData(imageCount).imageSegmenter = 'Raw Image';
                    imageData(imageCount).trained = 'NA';
                    imageData(imageCount).IOU = -2;
                    imageData(imageCount).TP = -2;
                    imageData(imageCount).FP = -2;
                    
                    imageData(imageCount).image = imread(filename);
                    imageData(imageCount).minArea = 0;
                end

                

            end
            waitbar(image/length(imageDir))
        end
        
        cd ..\
        cd labels
        imageDir = dir;
        imageDir(1).folder
        for image = 1:length(imageDir)
            waitbar(0, h,"Loading Ground Truths")

           
            if ~imageDir(image).isdir
                filename = imageDir(image).name;
                if ~contains(filename, "MATLAB")
                    imageCount = imageCount+1;
                    
                    imageData(imageCount).filename = filename;
                    imageName = split(imageDir(image).name, '.');
                    imageData(imageCount).imageName = imageName{1};
                    if contains(imageDir(image).folder, '\')
                        folder = split(imageDir(image).folder, '\');
                    elseif contains(imageDir(image).folder, '/')
                        folder = split(imageDir(image).folder, '/');
                    end
                    imageData(imageCount).imageSource = folder{end-1};
                    imageData(imageCount).imageSegmenter = 'Ground Truth';
                    imageData(imageCount).trained = 'NA';
                    imageData(imageCount).IOU = -2;
                    imageData(imageCount).IOU_CB = -2;
                    imageData(imageCount).IOU_indivOrganoids_CB = -2;
                    imageData(imageCount).TP = -2;
                    imageData(imageCount).FP = -2;
                    
                    currImage = imread(filename);
                    if length(size(currImage)) == 3
                        currImage = rgb2gray(currImage);
                    end
                    if contains(imageData(imageCount).imageSource, "EB")
                        currImage = currImage==255;
                    end
                  
                    currImage = imopen(logical(currImage),ones(5));
                
                    imageData(imageCount).image = currImage;
                    statsBW = regionprops(bwconncomp(currImage), 'Area');
                    imageData(imageCount).minArea = min([statsBW.Area]);
                end
                

            end
            waitbar(image/length(imageDir))
        end
        cd ..\
        cd ..\
    end
        
end


extractorPattern = "pp_";
labelerPattern = "";
segPattern = "_colored_1";
IDpattern = "_id-labeled";
outputPatterns = containers.Map(["OrganoSeg2", "OrganoSeg1", "OrganoLabeler","OrganoID","OrgaExtractor"],{segPattern, segPattern, labelerPattern, IDpattern, extractorPattern});
for result = 1:length(mainDir)
    
    if contains(mainDir(result).name, 'Results') & mainDir(result).isdir
        segmenter = split(mainDir(result).name, 'Results');
        segmenter = segmenter{1}; 
        
        cd(mainDir(result).name)   

        sourcesDir = dir;
        namePattern = outputPatterns(segmenter);
        for source = 1:length(sourcesDir)
            if contains(sourcesDir(source).name, 'Images')
                
                cd(sourcesDir(source).name)   
                sourceName = sourcesDir(source).name;                        
                
                if (strcmp(segmenter, "OrganoID") | strcmp(segmenter, "OrgaExtractor"))
                    cd trained
                end

                imageDir = dir;
                
                waitbar(0, h,sprintf("Loading %s, %s Results", segmenter, sourceName))

                for image = 1:length(imageDir)                    
                   
                    if ~imageDir(image).isdir & contains(imageDir(image).name, namePattern) & ~contains(imageDir(image).name, 'MATLAB')
                        imageCount = imageCount+1;
                        filename = imageDir(image).name;
                        imageData(imageCount).filename = filename;
                        imageName = split(imageDir(image).name, '.');
                        imageName = erase(imageName(1), namePattern);
                        imageData(imageCount).imageName = imageName{1};
                        
                        imageData(imageCount).imageSource = sourceName;
                        imageData(imageCount).imageSegmenter = segmenter;
                        if contains(segmenter, "OrganoSeg") | strcmp(segmenter, "OrganoLabeler")
                            imageData(imageCount).trained = 'NA';
                        elseif strcmp(segmenter, "OrganoID")
                            imageData(imageCount).trained = 'Yes';
                        elseif strcmp(segmenter, "OrgaExtractor")
                            imageData(imageCount).trained = 'Yes';
                        end
                        match = strcmp({imageData.imageName}, imageName{1}) & strcmp({imageData.imageSegmenter}, "Ground Truth");
                        segCPNumber = find(match);
                        segmentedCounterpartObject = imageData(segCPNumber);
                        segmentedCounterpart = segmentedCounterpartObject(1).image;
                        
                        currImage = imread(filename);
                        
                        

                        
                        if length(size(currImage)) == 3
                            currImage = rgb2gray(currImage);
                        end
                         

                        

                        if strcmp(segmenter, "OrganoLabeler")
                            currImage = imopen(currImage > 50, ones(5));
                        else
                            nonClearArea = imopen(logical(currImage),ones(5));
                            currImage(~nonClearArea) = 0;
                        end
                       
                        if ~contains(imageData(imageCount).imageSource, "OrganoLabeler")
                            nonClearArea = bwareaopen(logical(currImage), segmentedCounterpartObject(1).minArea);
                            currImage(~nonClearArea) = 0;
                        end
                        

                       
                       
                        imageData(imageCount).image = currImage;
                        statsBW = regionprops(bwconncomp(logical(currImage)), 'Area');
                        imageData(imageCount).minArea = min([statsBW.Area]);

                        
                        
                        
                        if ~strcmp(imageData(imageCount).imageSource, "OrganoLabelerImagesBrain")
                            imageData(imageCount).IOU_CB = getIOU(currImage, segmentedCounterpart, 1);
                            [IOU_array, TP, FP, FN, currCirc, currEcc, currSolid, currArea, circCP, eccCP, solidCP, areaCP] = getIndivIOU(currImage, segmentedCounterpart, 1, segmenter);
                            
                        else
                            imageData(imageCount).IOU_CB = getIOU(currImage, segmentedCounterpart, 0);
                            [IOU_array, TP, FP, FN, currCirc, currEcc, currSolid, currArea, circCP, eccCP, solidCP, areaCP] = getIndivIOU(currImage, segmentedCounterpart, 0, segmenter);
                        end

                        imageData(imageCount).IOU_indivOrganoids_CB = IOU_array;
                        imageData(imageCount).TP = TP;
                        imageData(imageCount).FP = FP;
                        
                        imageData(segCPNumber).Area = areaCP;



                        
                        

        
                    end
                    waitbar(image/length(imageDir))
                end
                if (strcmp(segmenter, "OrganoID") | strcmp(segmenter, "OrgaExtractor"))
                    cd ..\
                end
                if (strcmp(segmenter, "OrganoID") | strcmp(segmenter, "OrgaExtractor")) 
                                         
                    
                    cd default
                    imageDir = dir;
                    imageDir(1).folder
                    waitbar(0, h,sprintf("Loading %s - Default, %s Results", segmenter, sourceName))

                    for image = 1:length(imageDir)                    
                       
                        if ~imageDir(image).isdir & contains(imageDir(image).name, namePattern) & ~contains(imageDir(image).name, 'MATLAB')
                            imageCount = imageCount+1;
                            filename = imageDir(image).name;
                            imageData(imageCount).filename = filename;
                            imageName = split(imageDir(image).name, '.');
                            imageName = erase(imageName(1), namePattern);
                            imageData(imageCount).imageName = imageName{1};
                            
                            imageData(imageCount).imageSource = sourceName;
                            imageData(imageCount).imageSegmenter = segmenter;
                            imageData(imageCount).trained = 'No';
                            currImage = imread(filename);
                            match = strcmp({imageData.imageName}, imageName{1}) & strcmp({imageData.imageSegmenter}, "Ground Truth");
                            segCPNumber = find(match);
                            segmentedCounterpartObject = imageData(segCPNumber);
                            segmentedCounterpart = segmentedCounterpartObject(1).image;
                            

                            if length(size(currImage)) == 3
                                currImage = rgb2gray(currImage);
                            end
                            
                            nonClearArea = imopen(logical(currImage),ones(5));
                            currImage(~nonClearArea) = 0;
                            
                            
                           
                            if ~contains(imageData(imageCount).imageSource, "OrganoLabeler")
                                nonClearArea = bwareaopen(logical(currImage), segmentedCounterpartObject(1).minArea);
                                currImage(~nonClearArea) = 0;
                            end
                          

                            
                            imageData(imageCount).image = currImage;
                            
                            if ~strcmp(imageData(imageCount).imageSource, "OrganoLabelerImagesBrain")
                                imageData(imageCount).IOU_CB = getIOU(currImage, segmentedCounterpart, 1);
                                [IOU_array, TP, FP, FN, currCirc, currEcc, currSolid, currArea, circCP, eccCP, solidCP, areaCP] = getIndivIOU(currImage, segmentedCounterpart, 1, segmenter);
                                
                            else
                                imageData(imageCount).IOU_CB = getIOU(currImage, segmentedCounterpart, 0);
                                [IOU_array, TP, FP, FN, currCirc, currEcc, currSolid, currArea, circCP, eccCP, solidCP, areaCP] = getIndivIOU(currImage, segmentedCounterpart, 0, segmenter);
                            end
                            imageData(imageCount).IOU_indivOrganoids_CB = IOU_array;
                            imageData(imageCount).TP = TP;
                            imageData(imageCount).FP = FP;
                            
                            
                            imageData(segCPNumber).Area = areaCP;
                            


            
                        end
                        waitbar(image/length(imageDir))

                    end
                    cd ..\

                end
                cd ..\
            end
    
        end
        cd ..\
    end
            
end
close(h)

%% 

% match = strcmp({imageData.imageSource}, "OrganoLabelerImagesBrain") & strcmp({imageData.imageSegmenter}, "Ground Truth");
% organoLabelerMasks = imageData(find(match));
% for i = 1:length(organoLabelerMasks)
%     figure(100+i)
%     imshow(label2rgb(labelmatrix(bwconncomp(organoLabelerMasks(i).image)), 'hsv', 'k', 'shuffle')*0.9)
% end

%% 
sources = ["OrgaExtractorImages", "OrganoIDImagesLung", "OrganoIDImagesPDAC", "OrganoLabelerImagesBrain", "OrganoLabelerImagesEB", "OrganoSegImagesBreast"];
segmenters = ["OrgaExtractor", "OrganoID", "OrganoLabeler", "OrganoSeg1", "OrganoSeg2"];

%IOUData = cell(length(segmenters)+1, length(sources));
IOU_CBData = cell(length(segmenters)+1, length(sources));
IOU_IndivCBData = cell(length(segmenters)+1, length(sources));
TPData = cell(length(segmenters)+1, length(sources));
FPData = cell(length(segmenters)+1, length(sources));
labels = cell(length(segmenters)+1, length(sources));
sourceNum = 0;


for source = sources
    sourceNum = sourceNum + 1;
    segmenterNum = 0;
    subDataGT = getSubData(getSubData(imageData, "imageSource", source), "imageSegmenter", "Ground Truth");
    
    
    
    for segmenter = segmenters
        
        segmenterNum = segmenterNum + 1;
        
        if strcmp(segmenter, "OrganoID") | strcmp(segmenter, "OrgaExtractor")
            
            subDataTrained = getSubData(getSubData(getSubData(imageData, "trained", "Yes"), "imageSource", source), "imageSegmenter", segmenter);
            
            IOU_CBData{segmenterNum, sourceNum} = [subDataTrained(:).IOU_CB];
            TPData{segmenterNum, sourceNum} = [subDataTrained(:).TP];
            FPData{segmenterNum, sourceNum} = [subDataTrained(:).FP];
            IOU_Indiv = [];
            
            for subImage = 1:length(subDataTrained)
                IOU_Indiv = [IOU_Indiv subDataTrained(subImage).IOU_indivOrganoids_CB];
                
            end
            IOU_IndivCBData{segmenterNum, sourceNum} = IOU_Indiv;
            labels{segmenterNum, sourceNum} = strcat(segmenter, source);
            segmenterNum = segmenterNum + 1;
            

            if contains(source, segmenter) 
                subDataDefault = getSubData(getSubData(getSubData(imageData, "trained", "Yes"), "imageSource", source), "imageSegmenter", segmenter);
            else
                subDataDefault = getSubData(getSubData(getSubData(imageData, "trained", "No"), "imageSource", source), "imageSegmenter", segmenter);
            end
            
            IOU_CBData{segmenterNum, sourceNum} = [subDataDefault(:).IOU_CB];
            TPData{segmenterNum, sourceNum} = [subDataDefault(:).TP];
            FPData{segmenterNum, sourceNum} = [subDataDefault(:).FP];
            IOU_Indiv = [];
            
            for subImage = 1:length(subDataDefault)
                IOU_Indiv = [IOU_Indiv subDataDefault(subImage).IOU_indivOrganoids_CB];
                
            end
            IOU_IndivCBData{segmenterNum, sourceNum} = IOU_Indiv;
            labels{segmenterNum, sourceNum} = strcat(segmenter, source);
            
        else
            subData = getSubData(getSubData(imageData, "imageSource", source), "imageSegmenter", segmenter);
            
            IOU_CBData{segmenterNum, sourceNum} = [subData(:).IOU_CB];
            TPData{segmenterNum, sourceNum} = [subData(:).TP];
            FPData{segmenterNum, sourceNum} = [subData(:).FP];
            
            IOU_Indiv = [];
            
            for subImage = 1:length(subData)
                IOU_Indiv = [IOU_Indiv subData(subImage).IOU_indivOrganoids_CB];
                
            end
            IOU_IndivCBData{segmenterNum, sourceNum} = IOU_Indiv;
           
            labels{segmenterNum, sourceNum} = strcat(segmenter, source);
            
            
        end
    end
end




%%  KS TESTS
segmenterAndTrainedNames = ["OrgaExtractor - Trained", "OrgaExtractor - Default", "OrganoID - Trained","OrganoID - Default", "OrganoLabeler", "OrganoSeg1", "OrganoSeg2"];
comparisonNamesTable = cell(length(segmenterAndTrainedNames),length(segmenterAndTrainedNames));
numSegmented = zeros(length(segmenterAndTrainedNames), length(sources));

for source = 1:length(sources)
    adjSourcePTable = zeros(length(segmenterAndTrainedNames),length(segmenterAndTrainedNames));
    sourcePTable = zeros(length(segmenterAndTrainedNames),length(segmenterAndTrainedNames));
    sourceKSTable = zeros(length(segmenterAndTrainedNames),length(segmenterAndTrainedNames));
   
    
    for seg1=1:length(segmenterAndTrainedNames)
        
        numSegmented(seg1, source) = sum(~(isnan(IOU_IndivCBData{seg1, source})));

        for seg2 = 1:length(segmenterAndTrainedNames)
            
            comparisonNamesTable{seg1,seg2} = sprintf("%s vs. %s", segmenterAndTrainedNames(seg1), segmenterAndTrainedNames(seg2));
            [h,p,ksStat] = kstest2(IOU_IndivCBData{seg1, source}, IOU_IndivCBData{seg2, source});
            
            sourcePTable(seg1,seg2) = p;
            sourceKSTable(seg1,seg2) = ksStat;
            if contains(sources(source), "Extractor") | contains(sources(source), "ID")
                numComparisons = 5;
            else
                numComparisons = 6;
            end

            adjSourcePTable(seg1,seg2) = 1-(1-p)^numComparisons;
        end
    end
    adjPs.(sources(source)) = adjSourcePTable;
    ps.(sources(source)) = sourcePTable;
    ksStats.(sources(source)) = sourceKSTable;
   
end







%% 



% 
%% 




truePositive = [];
falsePositive = [];
 

namesForGroups = ["Colon - OrgaExtractor", "Lung - OrganoID", "PDAC - OrganoID", "Brain - OrganoLabeler", "Embryoid - OrganoLabeler", "Breast - OrganoSeg"];
for group = 1:length(namesForGroups)
    tempTP = zeros(1,length(segmenterAndTrainedNames));
    tempFP = zeros(1,length(segmenterAndTrainedNames));
    
    figure(group)
    segmenterData = zeros(length(IOU_IndivCBData{1, group}), length(segmenterAndTrainedNames));
    
        
    subDataGT = getSubData(getSubData(imageData, "imageSource", sources(group)), "imageSegmenter", "Ground Truth");
    
    area_Indiv = [];
    for subImage = 1:length(subDataGT)
        
        area_Indiv = [area_Indiv subDataGT(subImage).Area];
        
    end


    dataArray = [];
    ticksArray = cell(1,length(namesForGroups));
    segmenterVar = categorical(segmenterAndTrainedNames);
    for segmenter = 1:length(segmenterAndTrainedNames)
        
        segmenterData(:,segmenter) = [IOU_IndivCBData{segmenter, group}]';
        tempTP(segmenter) = sum(TPData{segmenter,group})/length(segmenterData);
        tempFP(segmenter) = sum(FPData{segmenter,group})/length(segmenterData);
        
        
    
        
    end
    truePositive = [truePositive;tempTP];
    falsePositive = [falsePositive;tempFP];
    area_Indiv = area_Indiv';
    
    
    y = repmat(segmenterVar,size(segmenterData, 1),1);
    if contains(sources(group), 'Brain')
        markerAlpha = 1;
    else
        markerAlpha = 0.3;
    end
    


    scatterplotobj = swarmchart(segmenterData,y, [], colors7, 'filled', 'o', 'MarkerFaceAlpha',markerAlpha, YJitter="density");
    %scatterplotobj = swarmchart(segmenterData,y, [], colors7, 'filled', 'o', 'MarkerFaceAlpha','flat', 'AlphaData', area_Indiv(:,1), YJitter="density");
    
    xlim([0, 1])
    
    %a = swarmchart(segmenterData,y, 50, colors7,  '.',  YJitter="density")

    hold on
    
    for segmenter = 1:length(segmenterAndTrainedNames)
        
       
        medIOU = median([IOU_IndivCBData{segmenter, group}], 'omitmissing');
        lowIOU = prctile([IOU_IndivCBData{segmenter, group}], 25);
        highIOU = prctile([IOU_IndivCBData{segmenter, group}], 75);

        numEntries = sum(~isnan([IOU_IndivCBData{segmenter, group}]));
        
        a = errorbar(medIOU, segmenterVar(segmenter), medIOU-lowIOU, highIOU-medIOU, 'horizontal', 'k|', 'MarkerSize', 12);
        hold on
    end
    pbaspect([1,1,1])
    ax = gca;
    hold off
    % exportFileName = sprintf("%s_indivOrg_jitter_final.pdf",sources(group))
    % exportgraphics(ax,exportFileName,'ContentType','vector')
    
   

    



    


end

%% Compare shared segmentations organoId organoSeg PDac
organoID_PDAC = IOU_IndivCBData{4, 3};
organoSeg2_PDAC = IOU_IndivCBData{7, 3};
segmenterVar = categorical(["OrganoID", "OrganoSeg2"]);
sharedInd = intersect(find(~isnan(organoID_PDAC)), find(~isnan(organoSeg2_PDAC)));
segmenterData = [organoID_PDAC(sharedInd)' organoSeg2_PDAC(sharedInd)'];
y = repmat(segmenterVar,size(segmenterData, 1),1);
figure(10)
scatterplotobj = swarmchart(segmenterData,y, [], [colors7(4,:);colors7(7,:)], 'filled', 'o', 'MarkerFaceAlpha',markerAlpha, YJitter="density");
hold on
errorbar(median(organoID_PDAC(sharedInd)), segmenterVar(1), median(organoID_PDAC(sharedInd)) - prctile(organoID_PDAC(sharedInd), 25), prctile(organoID_PDAC(sharedInd), 75)-median(organoID_PDAC(sharedInd)), 'horizontal', 'k|', 'MarkerSize', 12);
errorbar(median(organoSeg2_PDAC(sharedInd)), segmenterVar(2), median(organoSeg2_PDAC(sharedInd)) - prctile(organoSeg2_PDAC(sharedInd), 25), prctile(organoSeg2_PDAC(sharedInd), 75)-median(organoSeg2_PDAC(sharedInd)), 'horizontal', 'k|', 'MarkerSize', 12);
[h,p] = kstest2(organoID_PDAC(sharedInd), organoSeg2_PDAC(sharedInd));
hold off
pbaspect([1,1,1])
stylegraph(gca)

figure(11)
subDataGT = getSubData(getSubData(imageData, "imageSource", sources(3)), "imageSegmenter", "Ground Truth");

area_Indiv = [];
for subImage = 1:length(subDataGT)
    
    area_Indiv = [area_Indiv subDataGT(subImage).Area];
end
segmenterVar = categorical(["OrganoID", "OrganoSeg2", "Legend"]);

area_Indiv = area_Indiv';
area_Indiv = repmat(area_Indiv, 1,2);
legendCol = NaN(1,length(organoSeg2_PDAC));
legendCol(1:7) = [0.1 0.2 0.3 0.4 0.5 0.6 0.7];
legendColArea = zeros(1,length(organoSeg2_PDAC));
legendColArea(1:7) = [100 200 400 800 1600 3200 6400];
segmenterData = [organoID_PDAC' organoSeg2_PDAC' legendCol'];
area_Indiv = [area_Indiv legendColArea'];
y = repmat(segmenterVar,size(segmenterData, 1),1);


area_Indiv_inv = (area_Indiv.^-1);
scatterplotobj = swarmchart(segmenterData, y, area_Indiv_inv*20000, [colors7(4,:);colors7(7,:);[0.5 0.5 0.5]], 'filled', 'o', 'MarkerFaceAlpha',markerAlpha, 'MarkerEdgeAlpha', 0.5, 'MarkerEdgeColor', [0 0 0], YJitter="density");
pbaspect([1,1,1])
stylegraph(gca)
xlim([0 1])
%% 















images = [1,22,33, 72,73, 99];


dilateSizes = [10 9 13 3 3 11];

colorsSources = [48 213 200; 191 91 23; 191 91 23; 240 2 127; 240 2 127; 191 135 212]/255;
figure(20)
repImageTiles = tiledlayout(2,length(images));
for imageNum = 1:length(images)
    
    imageName = imageData(images(imageNum)).imageName;
    match = strcmp({imageData.imageName}, imageName) & strcmp({imageData.imageSegmenter}, "Raw Image");
    rawImage = imageData(find(match));
    
    
    fullImage = rawImage(1).image;
    fullImage = double(fullImage);
    fullImage = (fullImage - min(fullImage(:)))./(max(fullImage(:)) - min(fullImage(:)));
   

    match = strcmp({imageData.imageName}, imageName) & strcmp({imageData.imageSegmenter}, "Ground Truth");
    gtSeg = imageData(find(match));
    gtImage = gtSeg(1).image;
    
    minDim = min(size(fullImage(:,:,1)));
    source = split(imageData(images(imageNum)).imageSource, "Images");
    source = source{1};
     if contains(imageData(images(imageNum)).imageSource, "PDAC")
        fullImage = imlocalbrighten(fullImage, 0.8);
     elseif contains(source, "ID") 
         fullImage = imlocalbrighten(fullImage, 0.3);
     elseif contains(source, "Seg") 
         fullImage = imlocalbrighten(fullImage, 0.2);
     end
  

    match = strcmp({imageData.imageName}, imageName) & contains({imageData.imageSegmenter}, source) & ~strcmp({imageData.trained}, 'No');
    
    sourceImage = imageData(find(match));
    
    if strcmp(sourceImage(1).imageSegmenter, 'OrganoID')
        sourceImage = sourceImage(1).image;
        mag = imgradient(sourceImage);
        perim = bwskel(mag > 0);
        
    else
        sourceImage = sourceImage(1).image;
        perim = bwperim(sourceImage);       
        
    end
    dil_perim = imdilate(perim, ones(dilateSizes(imageNum)));
        

    match = strcmp({imageData.imageName}, imageName) & contains({imageData.imageSegmenter}, "Seg2");
    org2Image = imageData(find(match));
    org2Image = org2Image(1).image;

    if contains(source, 'Seg')
        fullImage = (fullImage(1:minDim,200:200+minDim));
        dil_perim = (dil_perim(1:minDim,200:200+minDim));
        org2Image = (org2Image(1:minDim,200:200+minDim));
        gtImage = (gtImage(1:minDim,200:200+minDim));
    else

        fullImage = (fullImage(1:minDim,1:minDim));
        dil_perim = (dil_perim(1:minDim,1:minDim));
        org2Image = (org2Image(1:minDim,1:minDim));
        gtImage = (gtImage(1:minDim,1:minDim));
    end
    
    imageWithGT = imoverlay(fullImage, imdilate(bwperim(gtImage), ones(dilateSizes(imageNum))), [0 0 0]);
   

    nexttile(imageNum)    
    imshow(imoverlay(imageWithGT, dil_perim, colorsSources(imageNum, :)))
    overlaySource = imoverlay(imageWithGT, dil_perim, colorsSources(imageNum, :));
    
    
    

    nexttile(imageNum+length(images))
    imshow(imoverlay(imageWithGT, imdilate(bwperim(org2Image), ones(dilateSizes(imageNum))), colors7(7,:)))
    overlayOrg2 = imoverlay(imageWithGT, imdilate(bwperim(org2Image), ones(dilateSizes(imageNum))), colors7(7,:));
    
end










%% FUNCTIONS



function IOU = getIOU(currImage, segmentedCounterpart, clearBorder)
    if clearBorder == 1
        currImage = logical(imclearborder(currImage));
        segmentedCounterpart = imclearborder(segmentedCounterpart);
    end
    try 
        intersection = sum(currImage & segmentedCounterpart);
        union = sum(currImage | segmentedCounterpart);
    catch
        currImage = imresize(currImage, size(segmentedCounterpart));
        intersection = sum(currImage & segmentedCounterpart);
        union = sum(currImage | segmentedCounterpart);
    end    
    IOU = intersection/union;
end

function showOverlay(imageData, imageName, num)


match = strcmp({imageData.imageName}, imageName) & strcmp({imageData.imageSegmenter}, "Ground Truth");
segmentedCounterpart = imageData(find(match));
segmentedCounterpart = segmentedCounterpart(1).image;
fullLabel = segmentedCounterpart;
match = strcmp({imageData.imageName}, imageName) & strcmp({imageData.imageSegmenter}, "Raw Image");
rawImage = imageData(find(match));
fullImage = rawImage(1).image;
fullImage = double(fullImage);
fullImage = (fullImage - min(fullImage(:)))./(max(fullImage(:)) - min(fullImage(:)));
match = strcmp({imageData.imageName}, imageName) & ~strcmp({imageData.imageSegmenter}, "Ground Truth") & ~strcmp({imageData.imageSegmenter}, "Raw Image");
allImages = imageData(find(match));

if length(allImages) == 5
    imageMap = {[0 0 1], [1 0 0], [1 0 1], [1 1 0], [0 1 0]};
else
    imageMap = {[0 0 1], [1 0 0], [1 1 0], [0 1 0]};
end
for image = 1:length(allImages)
    
    allImages(image).imageSegmenter
    
    try 
        perim = imdilate(bwperim(logical(allImages(image).image)),ones(3));
        if contains(allImages(image).imageSegmenter, "OrganoSeg") | contains(allImages(image).imageSegmenter, "Labeler")
            if contains(allImages(image).imageSegmenter, "OrganoSeg")
                organoSegIOU = allImages(image).IOU;
            else
                orgaLabelerIOU = allImages(image).IOU;
            end
        else
            perim = zeros(size(perim));
        end
        fullLabel = imoverlay(fullLabel, perim, imageMap{image});
        fullImage = imoverlay(fullImage, perim, imageMap{image});
    catch
        resized = imresize(logical(allImages(image).image), size(fullLabel(:,:,1)));
        perim = imdilate(bwperim(resized), ones(3));
        if contains(allImages(image).imageSegmenter, "OrganoSeg") | contains(allImages(image).imageSegmenter, "Labeler")
            if contains(allImages(image).imageSegmenter, "OrganoSeg")
                organoSegIOU = allImages(image).IOU;
            else
                orgaLabelerIOU = allImages(image).IOU;
            end
        else
            perim = zeros(size(perim));
        end
        fullLabel = imoverlay(fullLabel, perim, imageMap{image});
        fullImage = imoverlay(fullImage, perim, imageMap{image});
    end
    

end

nexttile(num)
imshow(fullLabel)
% hold on
% for i = 1:length(imageMap)
%     plot(0,0, 'Color', imageMap{i})
% end
% if length(allImages) == 4
%     
%     legend("OrgaExtractor", "OrganoID", "OrganoLabeler", "OrganoSeg")
% else
%     legend("OrgaExtractor", "OrganoID-Trained", "OrganoID-Default", "OrganoLabeler", "OrganoSeg")
% end
% hold off

nexttile(num+3)



imshow(fullImage, []);

title(sprintf("Seg %.2f Lab %.2f", organoSegIOU, orgaLabelerIOU))
% hold on
% for i = 1:length(imageMap)
%     plot(0,0, 'Color', imageMap{i})
% end
% if length(allImages) == 4
%     
%     legend("OrgaExtractor", "OrganoID", "OrganoLabeler", "OrganoSeg")
% else
%     legend("OrgaExtractor", "OrganoID-Trained", "OrganoID-Default", "OrganoLabeler", "OrganoSeg")
% end
% hold off

end
    
function subData = getSubData(data, variable, value)
match = strcmp({data.(variable)}, value);
subData = data(find(match));
end
    

function [IOU_array, TP, FP, FN, currCirc, currEcc, currSolid, currArea, circCP, eccCP, solidCP, areaCP] = getIndivIOU(currImage, segmentedCounterpart, clearBorder, segmenter)
    if clearBorder == 1
        currImage = imclearborder(currImage);
        segmentedCounterpart = imclearborder(segmentedCounterpart);
    end
    try 
        intersection = sum(currImage & segmentedCounterpart);        
    catch
        currImage = imresize(currImage, size(segmentedCounterpart));
    end
    if contains(segmenter, 'ID')
        currCC = getCCFromLabels(currImage);
    else
        currCC = bwconncomp(currImage);
    end
    segCC = bwconncomp(segmentedCounterpart);
    IOU_withSegIndices = zeros(1,length(segCC.PixelIdxList));
    matched_withCurrIndices = zeros(1,length(currCC.PixelIdxList));
    statsBWSegCurr = regionprops(currCC, 'Area', 'Circularity', 'Eccentricity', 'Solidity');
    circularity = [statsBWSegCurr.Circularity];
    eccentricity = [statsBWSegCurr.Eccentricity];
    solidity = [statsBWSegCurr.Solidity];
    area = [statsBWSegCurr.Area];
    currCirc = zeros(1,length(segCC.PixelIdxList));
    currEcc = zeros(1,length(segCC.PixelIdxList));
    currSolid = zeros(1,length(segCC.PixelIdxList));
    currArea = zeros(1,length(segCC.PixelIdxList));
    for segID = 1:length(segCC.PixelIdxList)
        tempIOU = nan;
        mostIntersect = 0;
        mostIntersectId = nan;
        segSpherPixels = segCC.PixelIdxList{segID};
        for currID = 1:length(currCC.PixelIdxList)
            orgIntersect = length(intersect(currCC.PixelIdxList{currID}, segSpherPixels));
            
            if orgIntersect > mostIntersect & orgIntersect > 0
                currIOU = orgIntersect/length(union(currCC.PixelIdxList{currID}, segSpherPixels));
                if matched_withCurrIndices(currID) == 0
                    mostIntersectId = currID;
                    mostIntersect = orgIntersect;
                    tempIOU = currIOU;
                else
                    otherSeg = matched_withCurrIndices(currID);
                    if currIOU > IOU_withSegIndices(otherSeg)
                        matched_withCurrIndices(currID) = 0;
                        IOU_withSegIndices(otherSeg) = nan;
                        mostIntersectId = currID;
                        mostIntersect = orgIntersect;
                        tempIOU = currIOU;                    
                    end
                end
                
                
            end
        end
        IOU_withSegIndices(segID) = tempIOU;
        if ~isnan(mostIntersectId)
            matched_withCurrIndices(mostIntersectId) = segID;
            currCirc(segID) = circularity(mostIntersectId);
            currEcc(segID) = eccentricity(mostIntersectId);
            currSolid(segID) = solidity(mostIntersectId);
            currArea(segID) = area(mostIntersectId);
        else
            currCirc(segID) = nan;
            currEcc(segID) = nan;
            currSolid(segID) = nan;
            currArea(segID) = nan;
        end
    end
    

    IOU_array = IOU_withSegIndices;
    
    totalManualSeg = length(segCC.PixelIdxList);
    totalAutoSeg = length(currCC.PixelIdxList);
    TP = sum(logical(matched_withCurrIndices));
    FP = length(matched_withCurrIndices) - sum(logical(matched_withCurrIndices));
    FN = length(IOU_array) - sum(logical(matched_withCurrIndices));
    statsBWSegCP = regionprops(segCC, 'Area', 'Circularity', 'Eccentricity', 'Solidity');
    
    circCP = [statsBWSegCP.Circularity];
    solidCP = [statsBWSegCP.Solidity];
    eccCP = [statsBWSegCP.Eccentricity];
    areaCP = [statsBWSegCP.Area];
end

function repImage = getRepImage(imageData, imageName, num, total)

    
    colors7 = [127 201 127; 48 213 200; 253 192 134; 255 255 40; 240,2,127; 191 175 212; 255, 194, 14]/255;

    match = strcmp({imageData.imageName}, imageName) & strcmp({imageData.imageSegmenter}, "Ground Truth");
    segmentedCounterpart = imageData(find(match));
    segmentedCounterpart = imclearborder(segmentedCounterpart(1).image);
    fullLabel = segmentedCounterpart;
    match = strcmp({imageData.imageName}, imageName) & strcmp({imageData.imageSegmenter}, "Raw Image");
    rawImage = imageData(find(match));
    fullImage = rawImage(1).image;
    fullImage = double(fullImage);
    fullImage = (fullImage - min(fullImage(:)))./(max(fullImage(:)) - min(fullImage(:)));
    fullImageOrig = fullImage;
    match = strcmp({imageData.imageName}, imageName) & ~strcmp({imageData.imageSegmenter}, "Ground Truth") & ~strcmp({imageData.imageSegmenter}, "Raw Image");
    allImages = imageData(find(match));
    source = allImages.imageSource
    [allImages.imageSegmenter]
    for machineLearner = ["OrganoID", "OrgaExtractor"]
        antimatch = contains({allImages.imageSource}, machineLearner) & contains({allImages.imageSegmenter}, machineLearner) & strcmp({allImages.trained}, "No")
        match = ~antimatch;
        allImages = allImages(find(match));
    end
    [allImages.imageSegmenter]
   
    
    sumMask = zeros(size(segmentedCounterpart));
    perim = nan;
    for image = 1:length(allImages)
    %     switch allImages(image).imageSegmenter
    %         case "OrgaExtractor"
    %             if strcmp(allImages(image).trained, 'Yes')
    %                 colorTemp = colors7(1,:);
    %             else
    %                 colorTemp = colors7(2,:);
    %             end
    %         case "OrganoID"
    % 
    %             if strcmp(allImages(image).trained, 'Yes')
    %                 colorTemp = colors7(3,:);
    %             else
    %                 colorTemp = colors7(4,:);
    %             end
    %         case "OrganoLabeler"
    %             colorTemp = colors7(5,:);
    %         otherwise
    %             colorTemp = colors7(6,:);
    %     end
    % 
    % 
    % 
    %     resizedTemp = imresize(allImages(image).image, size(fullLabel));
    %     perimTemp = imdilate(bwperim(resizedTemp), ones(3));
    %     fullImageTemp = imoverlay(fullImageOrig, perimTemp, [1 0 1]);
    % 
    %     fullImageTemp = fullImageTemp(:,1:size(fullImageTemp,1),:);
    %     size(fullImageTemp)
        % figure(2000)
        % imshow(fullImageTemp)
        % ax = gca;
        %exportgraphics(ax,sprintf("./figures/supplemental/%s_%s_%soverlay.pdf", allImages(image).imageSegmenter, allImages(1).imageSource,allImages(image).trained), 'ContentType', 'vector');
        
        try 
            mask = imclearborder(logical(allImages(image).image));
            if contains(allImages(image).imageSegmenter, "OrganoSeg")
                perim = imdilate(bwperim(allImages(image).image), ones(3));
            end
            sumMask = sumMask + mask;
            
        catch
            resized = imresize(allImages(image).image, size(fullLabel));
            mask = imclearborder(logical(resized));
            sumMask = sumMask + mask;
            
        end
        
            

    end
    sumMask(~segmentedCounterpart) = -sumMask(~segmentedCounterpart);
    
    
    max(sumMask(:))
    fullLabel = imoverlay(sumMask, perim, [0 0 1]);
    fullImage = imoverlay(fullImage, perim, [255/255 255/255 153/255]);
    %ax = nexttile(num);
    if length(allImages) == 7
        sumMask = sumMask/7;
        map = ([1 0 0; 6/7 0 0; 5/7 0 0; 4/7 0 0; 3/7 0 0;2/7 0 0;1/7 0 0; 0.5 0.5 0.5; 0 0 1/7; 0 0 2/7; 0 0 3/7; 0 0 4/7; 0 0 5/7;0 0 6/7;0 0 1; 255/255 255/255 153/255]);
        sumMask(perim) = 7/6;
        '7'
        allImages(image).imageName
    elseif length(allImages) == 6
        sumMask = sumMask/6;
         map = ([1 0 0; 5/6 0 0; 4/6 0 0; 3/6 0 0;2/6 0 0;1/6 0 0; 0.5 0.5 0.5; 0 0 1/6; 0 0 2/6; 0 0 3/6; 0 0 4/6; 0 0 5/6;0 0 1; 255/255 255/255 153/255]);
        sumMask(perim) = 1.2;
        '6'
        allImages(image).imageName
    end
    if size(sumMask,1)~=size(sumMask,2)
        if size(sumMask,1) > size(sumMask,2)
            sizeLim = size(sumMask,2);
        else
            sizeLim = size(sumMask,1);
        end
        sumMask = sumMask(1:sizeLim, 1:sizeLim);
        fullImage = fullImage(1:sizeLim, 1:sizeLim, :);
    end
    figure(100+num);
    imshow(sumMask)
    colormap(map)
    if length(allImages) == 6
        caxis([-1 7/6])
    else
        caxis([-1,1.2])
    end
%     figure(200+num)
%     imshow(allImages(end).image & segmentedCounterpart)
%     figure(300+num)
%     imshow(allImages(end).image & ~segmentedCounterpart)
%     figure(400+num)
%     imshow(~allImages(end).image & segmentedCounterpart)
%     figure(500+num)
%     imshow(imoverlay(fullImage, bwperim(segmentedCounterpart), [1 0 0]))

    ax = gca;
    %exportgraphics(ax,sprintf("./figures/finalCapstoneFigures_PostFilterSmall/%s_mask.pdf", allImages(1).imageSource ),'ContentType','vector')

    
    figure(200+num)
    imshow(fullImage, []);
    ax = gca;
    %exportgraphics(ax,sprintf("./figures/finalCapstoneFigures_PostFilterSmall/%s_overlay.pdf", allImages(1).imageSource ),'ContentType','vector')
    


end



function cc = getCCFromLabels(image)
    
    cc.Connectivity = 8;
    cc.ImageSize = size(image);
    
    pixelList = cell(1, max(image(:)));
    for i = 1:length(pixelList)
        pixelList{i} = find(image == i);
    end
    cc.PixelIdxList = pixelList;
    cc.PixelIdxList = cc.PixelIdxList(~cellfun('isempty',cc.PixelIdxList));
    cc.NumObjects = length(cc.PixelIdxList);

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