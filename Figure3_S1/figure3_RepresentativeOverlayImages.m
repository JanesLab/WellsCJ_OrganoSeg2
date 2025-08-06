if exist('perimCDFGoldInSegTile', 'var') == 0
    comparisonFigures;
end
dilateIndex = 1;
overlays = cell(1,length(hiddenProps.bulkI));
coords = cell(1,length(hiddenProps.bulkI));
centers = cell(1,length(hiddenProps.bulkI));
backgroundColor = [230 171 2]/255;

colors = {[55 126 184]/255, [77 175 74]/255, [255 192 203]/255, [255 127 0]/255, [166 86 40]/255, [152 78 163]/255};

%% 

for i = 1:length(hiddenProps.bulkI)
    image = hiddenProps.bulkI{i};
    image = double(image);
    image = (image - min(image(:)))./(max(image(:)) - min(image(:)));
    perimMask = zeros(size(image));
    perimMaskNoDil = zeros(size(image));
    goldCC = goldCCs{i};
    stats = regionprops(goldCC, "Centroid");
    centersOrig = [stats.Centroid];
    centers{i} = reshape(centersOrig, 2, length(stats))';
    
    
    for region = 1:length(goldCC.PixelIdxList)
        regionMask = zeros(size(image));
        regionMask(goldCC.PixelIdxList{region}) = 1;
        regionMask = bwperim(regionMask);
        regionMaskDil = imdilate(regionMask, ones(dilations(dilateIndex)));
        
        perimMask(regionMaskDil == 1) = 1;
        perimMaskNoDil(regionMask == 1) = 1;

    end
    
    


    
    imageWithGoldDil = imoverlay(image, perimMask, backgroundColor);
    imageWithGold = imoverlay(image, perimMaskNoDil, backgroundColor);

    coordsOrig = bwboundaries(logical(labelmatrix(ccs.borten{i})));
    coordsEdge = bwboundaries(logical(labelmatrix(ccs.edge{i})));
    
    badCoordsOrig = bwperim(logical(labelmatrix(ccs.borten{i}))) & ~perimMask;
    badCoordsEdge = bwperim(logical(labelmatrix(ccs.edge{i})))& ~perimMask;
    regCoordsOrig = bwperim(logical(labelmatrix(ccs.borten{i})));
    regCoordsEdge = bwperim(logical(labelmatrix(ccs.edge{i})));
    

    imageWithGoldDilOriginal = image;
    imageWithGoldDilEdge = image;
    imageWithGoldDilAdj = imageWithGoldDil;
    conditions = {image, imageWithGoldDilOriginal, imageWithGoldDilEdge, image};
    coords{i} = {0, coordsOrig, coordsEdge};
    regCoords{i} = {perimMask, regCoordsOrig, regCoordsEdge};
    badCoords{i} = {zeros(size(perimMask)), badCoordsOrig, badCoordsEdge};
    
    overlays{i} = conditions;
    

end
%% 
repImages = [4 5 7 8 10 16];
cropCoords = {[238,308,266,336], [178,298,362,482], [651, 751, 573, 673], [477,647,652,822],[522,622,372,472], [266,366,264,364]};
for imageInd = 1:length(repImages)
    imageNum = repImages(imageInd);
    regionID = find(mandersSegInGoldTotalPerim{imageNum}(1, :, dilateIndex) == mandersValues1(imageInd));
    
    
    
    conditionNames = {'Gold Standard Outline' 'Original Segmentation' 'Edge Corrected Segmentation'};
    
    goldCC = goldCCs{imageNum};
   
   
    
    
    

    figure(2+imageInd)
    clf
    perimCDFGoldInSegTile = tiledlayout(2,2);
    for condition = 1:3
        if condition == 1
            color = [166,206,227]/255;
        else
            color = colors{imageInd};
        end
        ax = nexttile;
        hold on
        if condition > 0
            currentCoords = coords{imageNum}{condition};

            if condition == 3
                lineStyle = ':';
            else
                lineStyle = '-';
            end
            currentCoords;
            
            lowPercentile = prctile(overlays{imageNum}{condition}(:),1);
            highPercentile = prctile(overlays{imageNum}{condition}(:),99)*1.05;
            if highPercentile>1
                highPercentile = 1;
            end

            if condition ~= 1
                overlay = imoverlay(imadjust(overlays{imageNum}{condition}, [lowPercentile;highPercentile] ), imdilate(regCoords{imageNum}{condition},ones(3)), color);
            else
                overlay = imoverlay(imadjust(overlays{imageNum}{condition}, [lowPercentile;highPercentile]), regCoords{imageNum}{condition}, color);
            end
            crop = cropCoords{imageInd};
           
            overlay = overlay(crop(1):crop(2), crop(3):crop(4), :);

            imshow(overlay)
             
            hold on
            
        end    
        
        hold off
       
    end
    title(perimCDFGoldInSegTile, strcat(['Image #', num2str(imageNum),':', ' ', 'Dilation = ',' ', num2str(dilations(dilateIndex)),' ', 'Pixels']))
end


