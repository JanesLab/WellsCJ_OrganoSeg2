getSegmentedData;
combineGoldStandard;
numImages = length(goldCCs);


mandersSegInGoldTotalPerim = cell(1, numImages);


goldCirc = cell(1,numImages);
segCircTotal = cell(1,numImages);
goldSol = cell(1,numImages);
segSolTotal = cell(1,numImages);
dilations = [5];
goldStandardPerims = cell(1,numImages);
example = 1;
for i = 1:length(goldCCs)
    goldCC = goldCCs{i};
    goldStats = regionprops(goldCC, 'Circularity', 'Solidity');
    goldCirc{i} = [goldStats.Circularity];
    goldSol{i} = [goldStats.Solidity];
    bortenCC = ccs.borten{i};
    edgSol = ccs.edge{i};
  
    currCCs = {bortenCC edgSol};

    
    segCircCurr = zeros(4, length(goldCC.PixelIdxList));
    segSolCurr = zeros(4, length(goldCC.PixelIdxList));
    

    numDilationLevels = length(dilations);
    mandersSegInGoldCurrPerim = zeros(4, length(goldCC.PixelIdxList), numDilationLevels);
   


    for cc = 1:length(currCCs)
        currCC = currCCs{cc};
        currCC.PixelIdxList = currCC.PixelIdxList(~cellfun('isempty',currCC.PixelIdxList));
        currCC.NumObjects = length(currCC.PixelIdxList);
        
        
        currStats = regionprops(currCC, 'Solidity', 'Circularity');
        
        currCirc = [currStats.Circularity];
        currSol = [currStats.Solidity];
        
        
        for goldSpher = 1:length(goldCC.PixelIdxList)
            
            mostIntersect = 0;
            mostIntersectId = nan;
            goldSpherPixels = goldCC.PixelIdxList{goldSpher};
            for currSpher = 1:length(currCC.PixelIdxList)
                spherIntersect = length(intersect(currCC.PixelIdxList{currSpher}, goldSpherPixels));
                if spherIntersect > mostIntersect
                    [M1, M2] = manders(currCC.PixelIdxList{currSpher}, goldSpherPixels);
                    if M1 > 0.5 && M2 > 0.5
                        mostIntersectId = currSpher;
                        mostIntesect = spherIntersect;
                    end
                end
            end

            if ~isnan(mostIntersectId)

                segSpherPixels = currCC.PixelIdxList{mostIntersectId};
    
    
                
                

                


                segCircCurr(cc, goldSpher) = currCirc(mostIntersectId);
                segSolCurr(cc, goldSpher) = currSol(mostIntersectId);


                
                for dilateNumber  = 1:numDilationLevels
                    dilateFactor = dilations(dilateNumber);
                    segSpherPixelsPerim = getPerimPixels(segSpherPixels, currCC.ImageSize, 0);
                    goldSpherPixelsPerim = getPerimPixels(goldSpherPixels, goldCC.ImageSize, dilateFactor);
                    [M1Perim, M2Perim] = manders(segSpherPixelsPerim, goldSpherPixelsPerim);
    
    
                    mandersSegInGoldCurrPerim(cc, goldSpher, dilateNumber) = M1Perim;
                end

                


            else
                
                segCircCurr(cc, goldSpher) = nan;
                segSolCurr(cc, goldSpher) = nan;

                mandersSegInGoldCurrPerim(cc, goldSpher, :) = nan;
               

                
            end


        end

        

    end
    
    segCircTotal{i} = segCircCurr;
    segSolTotal{i} = segSolCurr;

    
    mandersSegInGoldTotalPerim{i} = mandersSegInGoldCurrPerim;
    
       

end

%% 



mandersSegInGoldAvgPerim = zeros(4, length(mandersSegInGoldTotalPerim), numDilationLevels);

for i = 1:length(mandersSegInGoldTotalPerim)
    
   
    mandersSegInGoldAvgPerim(:, i, :) = mean(mandersSegInGoldTotalPerim{i},2,"omitnan");
end

figure(1)
tiledlayout(1,2)
colors = [191 135 212; 255, 194, 14]/255;
combinedGoldCirc = cell2mat(goldCirc);
combinedSegCirc = cell2mat(segCircTotal);
combinedGoldSol = cell2mat(goldSol);
combinedSegSol = cell2mat(segSolTotal);
circCorr = zeros(1,2);
solCorr = zeros(1,2);
for i = 1:2
    nexttile(1)
    circularityityCurr = combinedSegCirc(i,:); 
    solidityCurr = combinedSegSol(i,:); 
    notNanInd = find(~isnan(circularityityCurr));
    scatter(combinedGoldCirc(notNanInd), circularityityCurr(notNanInd), 20,  colors(i,:), 'o','filled');
    hold on
    [corr, p] = corrcoef(combinedGoldCirc(notNanInd), circularityityCurr(notNanInd));
    xlim([0 1])
    ylim([0 1])
    circCorr(i) = corr(1,2);
    
    pbaspect([1,1,1])
    stylegraph(gca)
    title("Circularity")
    
    nexttile(2)
    scatter(combinedGoldSol(notNanInd), solidityCurr(notNanInd), 20,  colors(i,:), 'o','filled');
    hold on
    title("Solidity")
    [corr, p] = corrcoef(combinedGoldSol(notNanInd), solidityCurr(notNanInd));
    xlim([0.5 1])
    ylim([0.5 1])
    solCorr(i) = corr(1,2);
    pbaspect([1,1,1])
    stylegraph(gca)


end
nexttile(1)
legend(sprintf("OrganoSeg: r = %.2f", circCorr(1)), sprintf("OrganoSeg2: r = %.2f", circCorr(2)))
hold off
nexttile(2)
legend(sprintf("OrganoSeg: r = %.2f", solCorr(1)), sprintf("OrganoSeg2: r = %.2f", solCorr(2)))
hold off


function [M1, M2] = manders(segment, gold)
    M1 = length(intersect(segment, gold))/length(segment);
    M2 = length(intersect(segment, gold))/length(gold);

end

function perimPixels = getPerimPixels(pixels, imageSize, dilateFactor)


mask = zeros(imageSize);
mask(pixels) = 1;
perim = bwperim(mask);
if dilateFactor ~= 0
    perim = imdilate(perim, ones(dilateFactor));
end
perimPixels = find(perim);
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
