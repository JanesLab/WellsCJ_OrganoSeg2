%%Used to generate trackingAnalysisResults.mat by applying different registration methods follwed by organoid tracking
%recommended to run on high-powered computing


D03_Center93 = [24 24 29 35 31 33; % 
    97 107 117 128 126 121;
    42 49 55 59 58 60;
    74 82 90 98 97 98
    49 56 64 66 66 68;
    61 68 77 80 80 81;
    nan nan 72 76 nan 77;
    39 40 47 52 50 51;
    94 101 110 121 118 119;
    56 60 68 71 72 73;
    34 nan 42 49 46 47;
    41 48 53 57 53 55];


D03_Edge93 = [86 94 102 112 111 109; % 
    98 113 124 136 134 139;
    112 128 141 156 153 160;
    118 136 148 164 162 170;
    119 137 149 166 163 171;
    126 142 157 172 168 176;
    nan 148 163 180 178 188;
    nan 149 164 177 175 183;
    129 150 165 181 179 189;
    101 116 128 143 142 146;
    113 129 142 157 154 161;
    nan 109 120 130 128 129];

B02_Center92 = [145 176 178 192 187 173;  
    223 274 290 305 297 269;
    310 370 386 410 396 368;
    269 320 321 343 339 309;
    137 156 168 185 182 170;
    322 389 404 421 402 378;
    207 267 272 287 291 254;
    306 355 358 381 374 345;
    142 171 174 191 186  178;
    260 313 323 330 317 312;
    268 323 325 345 343 319;
    236 291 295 322 313 286
    ];

B02_Edge92 = [350 424 441 460 429 386;  
    222 273 273 291 268 229;
    313 379 392 413 400 373;
    289 353 352 374 363 333;
    243 290 289 313 304 nan;
    278 331 341 368 359 332;
    344 421 442 470 456 431;
    319 388 405 431 420 393;
    234 285 281 309 298 267;
    nan 284 282 303 nan nan;
    242 293 293 317 300 263;
    241 287 298 321 314 278;
    227 278 278 310 303 272];



B03_Center87 = [5 18 22 nan 14 nan;   
    8 22 27 23 16 10;
    12 29 32 28 21 12;
    13 26 29 32 24 13;
    7 20 26 24 18 nan;
    6 19 21 18 12 7;
    14 32 39 37 30 22;
    nan nan 34 34 26 17];

B03_Edge87 = [
    nan 2 3 2 2 nan;
    nan 11 13 11 8 5;
    nan 5 6 5 5 2;
    nan 6 7 6 nan nan;
    nan 4 5 4 4 nan;
    nan 3 4 3 3 nan;
    nan nan nan 13 9 4];


B02_Edge87 = [11 25 19 23 25 19;   
    14 28 22 26 27 nan;
    12 26 nan nan nan nan;
    13 27 nan nan nan nan;
    16 31 24 28 28 20;
    8 21 16 nan 23 15
    15 29 21 25 nan nan];

B02_Center87 = [3 4 nan 3 nan nan;  
    nan 2 nan 5 4 2;
    nan nan 2 4 3 nan;
    nan 8 8 10 9 6;
    4 6 4 7 6 4;
    nan nan 6 9 7 5;
    nan 14 nan nan 13 8;
    nan 17 11 14 14 9;
    nan nan nan 17 20 14;
    nan 20 14 16 18 11];


B02_Center85 = [188 196 186 169 153 121;
    299 301 290 259 231 209;
    345 359 348 311 292 259;
 
 
    182 166 163 150 136 126;
    296 300 287 251 229 208;
    408 395 376 338 312 274;
    237 247 242 213 191 nan;
    252 251 241 211 194 173;
    353 350 335 298 282 247;
    183 187 180 160 145 130;
    195 204 191 177 162 152;
    296 300 287 251 229 208;
    295 297 288 252 228 206];

B02_Edge85 = [477 458 433 388 364 nan;
    491 467 451 409 377 335;
    492 470 455 412 375 336;
    469 450 430 389 365 324;
   
    484 465 444 401 369 nan;
    456 438 419 377 348 309;
    452 435 416 373 349 311;
    376 372 355 318 298 257;
    454 434 418 372 353 315;
    494 471 453 411 378 338;
    478 463 438 398 368 326;
    485 466 448 404 372 332];



C04_Center81 = [34 52 54 48 45 38;   
    82 156 172 161 162 135;
    nan 151 173 152 158 132;
    106 199 216 203 199 158;
    108 186 212 197 189 nan;
    121 217 224 207 202 161;
    155 272 298 278 283 223
    nan nan 292 265 268 nan;
    190 316 333 314 318 243;
    62 120 130 117 114 106;
    nan 162 169 154 150 nan;
    nan 148 152 138 133 nan;
    33 67 69 59 54 36;
    nan 96 96 86 78 61;
    42 87 85 72 64 50];
 
C04_Edge81 = [130 225 238 223 222 178;   
    148 261 284 260 269 206;
    209 344 374 351 360 275;
    176 299 320 299 306 nan;
    166 289 314 291 289 nan;
    183 303 324 297 292 213;
    195 320 347 329 331 256;
    nan 354 382 358 367 281;
    174 296 310 296 290 235;
    nan 206 211 205 200 160;
    206 336 335 315 nan nan;
    185 306 330 304 313 nan;
    140 250 266 239 241 184;
    193 321 344 322 327 253];

D02_Edge118 = [87 89 73 36 34 27;
    104 105 84 44 43 33;
   
    158 141 124 72 68 58;
    166 176 144 nan nan nan;
    nan 126 112 61 nan 52;
    nan 149 129 73 72 62;
    129 125 110 59 61 46;    
    32 50 33 20 18 18];

D02_Center118 = [nan 179 141 82 81 69;   
    nan 165 139 80 79 67;
    168 183 152 93 91 77;
    137 136 122 68 65 55;
    180 199 159 95 89 nan;
    39 46 26 15 16 13;
    nan 216 179 104 94 74;
    63 75 59 32 29 23;
    45 55 48 27 24 21;
    nan 40 35 16 14 16;
    17 16 12 5 4 6;
    110 109 96 50 55 47;
    nan nan 99 49 53 45];

B02_Edge118 = [61 65 84 71 79 65;  
    67 79 92 77 85 70;
    nan 60 79 66 75 61;
    44 45 59 51 64 54;
    58 64 82 70 78 63;
    64 72 88 76 84 69;
    nan 77 93 78 86 71;
    52 56 73 63 71 nan;
    71 nan 99 83 93 77;
    33 37 52 48 55 46];

B02_Center118 = [43 41 56 50 58 48;  
    nan 42 57 52 59 50;
    nan 34 46 43 46 40;
    nan 28 40 38 40 36;
    25 22 34 34 34 30;
    60 66 83 72 80 64;
    48 52 66 57 67 55;
    nan 54 68 59 68 53;
    50 53 67 58 69 56;
    57 nan nan 69 77 62;
    nan 15 27 29 31 28];

E02_Center116 = [33 33 41 36 33 28;  
    28 28 37 37 37 nan; 
    47 51 58 56 54 50;
    80 90 98 103 102 101;
    87 97 103 107 104 109;
    29 32 39 nan nan nan;
    45 48 54 54 51 53;
    52 63 64 62 61 61;
    70 76 80 70 69 71];

 
E02_Edge116 = [42 47 51 51 49 54;  
    58 66 68 68 66 67;
    95 98 104 89 86 90;
    79 87 92 95 93 94;
    76 78 83 90 88 88;
    67 75 81 83 80 79;
    49 53 60 59 55 58;
    30 30 40 38 35 33;
    17 14 18 nan nan nan;
    13 12 14 13 nan nan];

B02_Center116 = [103 108 104 107 96 89;   
    35 44 47 51 47 41;    
    46 56 52 56 55 51;
    44 52 54 57 56 52;
    89 100 96 98 90 81;
    71 81 81 84 79 73;
    85 96 93 96 87 82;
    27 33 37 40 43 40;
    38 46 48 53 nan nan;
    93 103 99 101 91 86];

B02_Edge116 = [
    117 125 123 126 117 nan;
    108 115 106 112 103 94;
    81 91 89 91 84 79;
    nan 77 77 83 78 71;
    76 87 88 nan 80 76;
    82 89 85 88 83 77;
    118 126 124 127 118 105;
    128 140 138 138 129 113;
    131 143 142 140 131 nan;
    97 106 102 108 97 nan;
    115 123 121 125 116 103;
    111 117 114 120 110 99];



%% 

manualTrackings = {B02_Edge116, B02_Center116; E02_Edge116, E02_Center116; B02_Edge118, B02_Center118; D02_Edge118, D02_Center118; C04_Edge81, C04_Center81; B02_Edge85, B02_Center85; B02_Edge87, B02_Center87; B03_Edge87,B03_Center87; B02_Edge92, B02_Center92; D03_Edge93, D03_Center93};

cases = [116, 116, 118, 118, 81, 85, 87, 87, 92, 93];
wellNums = ["B02", "E02", "B02", "D02", "C04", "B02", "B02", "B03","B02","D03"];
registrations = ["Rigid","Rigid-Rotate", "Non-rigid"];
smoothingDegrees = [0.5,1,2,3];
imageMethods = ["Centroids", "Image", "Mask"];


automatedTrackins = cell(1,length(wellNums));


filenames = {dir("/project/g_bme-janeslab/Najwa/Organoid Images - EVOS/Tamoxifen Cases - Tracking").name};
filenames
nextIndex = 0;
count = 0;
loadedProps = cell(1, length(wellNums));
correctComparisons = cell(1, length(wellNums));
totalComparisons = cell(1, length(wellNums));
percentRetainedNonRig = cell(1, length(wellNums));

percentCorrectRig = zeros(1, length(wellNums));
percentRetainedRig = cell(1, length(wellNums));
adjustedCentersAllSetsRig = cell(1, length(wellNums));
matchingIDsAllSetsRig = cell(1, length(wellNums));
totalTformAllSets = cell(1, length(wellNums));


adjustedCentersAllSetsNonRig = cell(1, length(wellNums));
matchingIDsAllSetsNonRig = cell(1, length(wellNums));
displacementFieldsAllSets = cell(1, length(wellNums));

centersAllSets = cell(1, length(wellNums));
bw_perim_dil_AllSets = cell(1, length(wellNums));
trackingResults = cell(1,length(wellNums));
parfor imageSet = 1:length(wellNums)
    fileInd = contains(filenames, wellNums(imageSet)) & contains(filenames, "BCO"+cases(imageSet)) & contains(filenames, "Tracking");
    loadedProps{imageSet} = load("/project/g_bme-janeslab/Najwa/Organoid Images - EVOS/Tamoxifen Cases - Tracking/" + filenames{fileInd}, 'hiddenProps');
    hiddenProps = loadedProps{imageSet}.hiddenProps;
    wellNums(imageSet) 
    
    ccsAll = cell(1,length(hiddenProps.bulkI));
    centersAll = cell(1,length(hiddenProps.bulkI));
    bw_perim_dil_All = cell(1,length(hiddenProps.bulkI));
    

    for i = 1:length(hiddenProps.bulkI)
        
        CClabel = labelmatrix(hiddenProps.cc{i});
                

        bwMock = logical(CClabel);
        bw_holesTemp = bwMock;
       
            
        %Find distance transform and smooth
        dist = bwdist(imcomplement(double(bw_holesTemp)));
        
        
        smooth = imgaussfilt(dist,5);
        
        %smooth = dist;
        
        dist2 = -smooth;
        dist2(~bw_holesTemp) = Inf;

        %Watershed based on distance transform, splitting organoids
        L = watershed(dist2);
        bw_holesTemp(L==0) = 0;
        %Remove noise based on size threshold
        bwlil = bwareaopen(bw_holesTemp,150);
        %Clear border
        bw_holesTemp = imclearborder(bwlil);
        bw_holesTemp = logical(bw_holesTemp);
        
        bwMock = bw_holesTemp;
        %bwMock(bwMock > 0) = 1;
        
        
        

        %Find perimeter
        bw_perim = double(bwperim(bwMock));
        %Make visible (dilate)
        bw_perim_dil_All{i} = imdilate(bw_perim,ones(3,3));

        cc = bwconncomp(bwMock);
        ccsAll{i} = cc;
        statsBW = regionprops(cc, 'Centroid');
        centersOrig = [statsBW.Centroid];
        centers = reshape(centersOrig, 2, length(statsBW))';
        centersAll{i} = centers;
       

        
        %%%%%%%%Edge
    end
    tracesEdge = manualTrackings{imageSet, 1};
    tracesCentral = manualTrackings{imageSet, 2};    
    combinedTraces = [tracesEdge; tracesCentral];
    thisCaseResults = struct;
    thisCaseResults.cc = ccsAll;
    thisCaseResults.centers = centersAll;
    thisCaseResults.combinedTraces = combinedTraces;
    thisCaseResults.wellNum = wellNums(imageSet);
    thisCaseResults.caseNum = cases(imageSet);
    
    for imageMethod = imageMethods
        imageMethod
        thisMethodResults = struct;
        for registration = registrations
            registration
            if strcmp(registration, "Non-rigid")
                
                for smoothing = smoothingDegrees
                    smoothing
                    thisInstanceResults = struct;
                    [adjustedCenters, transforms] = nonRigidRegister(ccsAll, centersAll, hiddenProps.bulkI, imageMethod, 100, smoothing);
                    matchingIDs = trackOrganoids(adjustedCenters, 50);
                    [correctComparisons, totalComparisons, correctRetained, totalRetained] = evaluatePerformance(matchingIDs, combinedTraces);
                    if smoothing == 0.5
                        totalName = "NonRigid_Smooth05"
                    else
                        totalName = "NonRigid_Smooth" + num2str(smoothing);
                    end
                    thisInstanceResults.adjustedCenters = adjustedCenters;
                   % thisInstanceResults.transforms = transforms;
                    thisInstanceResults.matchingIDs = matchingIDs;
                    thisInstanceResults.correctComparisons = correctComparisons;
                    thisInstanceResults.correctRetained = correctRetained;
                    thisInstanceResults.totalRetained = totalRetained;
                    thisInstanceResults.totalComparisons = totalComparisons;
                    thisMethodResults.(totalName) = thisInstanceResults;
                end
            else
               if contains(registration, "Rotate")
                    [adjustedCenters, transforms] = rigidRegister(ccsAll, centersAll, hiddenProps.bulkI, imageMethod, true);
                    
                    totalName = "RigidRotate";
                else
                    [adjustedCenters, transforms] = rigidRegister(ccsAll, centersAll, hiddenProps.bulkI, imageMethod, false);
                    
                    totalName = "Rigid";
                end
                matchingIDs = trackOrganoids(adjustedCenters, 50);
                [correctComparisons, totalComparisons, correctRetained, totalRetained] = evaluatePerformance(matchingIDs, combinedTraces);
                
                thisInstanceResults = struct;
                thisInstanceResults.adjustedCenters = adjustedCenters;
               % thisInstanceResults.transforms = transforms;
                thisInstanceResults.matchingIDs = matchingIDs;
                thisInstanceResults.correctComparisons = correctComparisons;
                thisInstanceResults.correctRetained = correctRetained;
                thisInstanceResults.totalRetained = totalRetained;
                thisInstanceResults.totalComparisons = totalComparisons;
                thisMethodResults.(totalName) = thisInstanceResults;
            end  
            
        end
        thisCaseResults.(imageMethod) = thisMethodResults;
    end
    trackingResults{imageSet} = thisCaseResults;
end



save("trackingAnalysisResultsWithRotate.mat", "trackingResults")





%% 




cases = [116, 116, 118, 118, 81, 85, 87, 87, 92, 93];
wellNums = ["B02", "E02", "B02", "D02", "C04", "B02", "B02", "B03","B02","D03"];
registrations = ["Rigid","Rigid-Rotate", "Non-rigid"];
smoothingDegrees = [0.5,1,2,3];
imageMethods = ["Centroids", "Image", "Mask"];




%trackingResultsTemp = trackingResults;
filenames = {dir("/project/g_bme-janeslab/Najwa/Organoid Images - EVOS/Tamoxifen Cases - Tracking").name};
filenames
nextIndex = 0;
count = 0;
loadedProps = cell(1, length(wellNums));
correctComparisons = cell(1, length(wellNums));
totalComparisons = cell(1, length(wellNums));
percentRetainedNonRig = cell(1, length(wellNums));

percentCorrectRig = zeros(1, length(wellNums));
percentRetainedRig = cell(1, length(wellNums));
adjustedCentersAllSetsRig = cell(1, length(wellNums));
matchingIDsAllSetsRig = cell(1, length(wellNums));
totalTformAllSets = cell(1, length(wellNums));


adjustedCentersAllSetsNonRig = cell(1, length(wellNums));
matchingIDsAllSetsNonRig = cell(1, length(wellNums));
displacementFieldsAllSets = cell(1, length(wellNums));

centersAllSets = cell(1, length(wellNums));
bw_perim_dil_AllSets = cell(1, length(wellNums));


for imageSet = 1:length(wellNums)
    
    

    
   
    
    for imageMethod = imageMethods
        
        thisMethodResults = struct;
        for registration = registrations
            
            if strcmp(registration, "Non-rigid")
                
                for smoothing = smoothingDegrees
                     if smoothing == 0.5
                        totalName = "NonRigid_Smooth05"
                    else
                        totalName = "NonRigid_Smooth" + num2str(smoothing);
                     end
                     
                    
                    
                    
                    [correctComparisons, totalComparisons, correctRetained, totalRetained] = evaluatePerformance(trackingResultsTemp{imageSet}.(imageMethod).(totalName).matchingIDs, trackingResultsTemp{imageSet}.combinedTraces);
                   
                    
                    trackingResultsTemp{imageSet}.(imageMethod).(totalName).correctComparisons = correctComparisons;
                    trackingResultsTemp{imageSet}.(imageMethod).(totalName).correctRetained = correctRetained;
                    trackingResultsTemp{imageSet}.(imageMethod).(totalName).totalRetained = totalRetained;
                    trackingResultsTemp{imageSet}.(imageMethod).(totalName).totalComparisons = totalComparisons;
                    
                end
            else
               if contains(registration, "Rotate")
                    
                    totalName = "RigidRotate";
                else
                    
                    
                    totalName = "Rigid";
                end
                
                [correctComparisons, totalComparisons, correctRetained, totalRetained] = evaluatePerformance(trackingResultsTemp{imageSet}.(imageMethod).(totalName).matchingIDs, trackingResultsTemp{imageSet}.combinedTraces);
                
               
                trackingResultsTemp{imageSet}.(imageMethod).(totalName).correctComparisons = correctComparisons;
                trackingResultsTemp{imageSet}.(imageMethod).(totalName).correctRetained = correctRetained;
                trackingResultsTemp{imageSet}.(imageMethod).(totalName).totalRetained = totalRetained;
                trackingResultsTemp{imageSet}.(imageMethod).(totalName).totalComparisons = totalComparisons;
            end  
            
        end
        
    end
end



save("trackingAnalysisResultsWithRotate_New.mat", "trackingResultsTemp")
%% 
evaluatePerformance(trackingResultsTemp{1}.("Mask").("Rigid").matchingIDs, trackingResultsTemp{1}.combinedTraces)







%% Functions


function [adjustedCenters, displacementFields] = nonRigidRegister(cc, centers, bulkI, regMethod, numIterations, smoothing)
    adjustedCenters = cell(1,length(bulkI));
    displacementFields = cell(1,length(bulkI));
    registeredBulkI = cell(1,length(bulkI));
    adjustedCenters{1} = centers{1};
    firstImageCCLabel = labelmatrix(cc{1});
    
    if strcmp(regMethod, "Mask")
        registeredBulkI{1} = logical(firstImageCCLabel);
        
        
    elseif strcmp(regMethod, "Image")
   
        registeredBulkI{1} = bulkI{1};
        %app.RegisteredBulkI{app.1} = imbinarize(app.BulkI{app.1}, 'adaptive', 'ForegroundPolarity', 'Dark', 'Sensitivity', 0.6);
    
        
    else
        tempBulkI = zeros(size(firstImageCCLabel));
        for center = 1:size(centers{1}, 1)
            tempBulkI(round(centers{1}(center, 2)), round(centers{1}(center, 1))) = 1;
        end
        registeredBulkI{1} = imdilate(tempBulkI, ones(3));
        

    end
    displacementFields{1} = zeros([size(bulkI{1}), 2]);
    firstPerim = bwperim(firstImageCCLabel);
    firstPerim = imdilate(firstPerim, ones(3,3));
    registeredPerimeter{1} = firstPerim;
    
    [meshX, meshY] = meshgrid(1:size(firstImageCCLabel,2),1:size(firstImageCCLabel,1));
    
    
    for selected = 2:length(bulkI)
        prev = selected - 1;
    %                 
    
        %firstImageCCLabel = labelmatrix(app.CC{prev});
        secondImageCCLabel = labelmatrix(cc{selected});
        
        %Perform registration with the previously registered image
        %so that displacement fields for all images align
        
        baseImage = registeredBulkI{prev};
        
        

        if strcmp(regMethod, "Mask")
            currImage = logical(secondImageCCLabel);
            
            
        elseif strcmp(regMethod, "Image")
       
            currImage = bulkI{selected};
            %app.RegisteredBulkI{app.1} = imbinarize(app.BulkI{app.1}, 'adaptive', 'ForegroundPolarity', 'Dark', 'Sensitivity', 0.6);
        
            
        else
            tempBulkI = zeros(size(secondImageCCLabel));
            for center = 1:size(centers{selected}, 1)
                tempBulkI(round(centers{selected}(center, 2)), round(centers{selected}(center, 1))) = 1;
            end
            currImage = imdilate(tempBulkI, ones(5));
            
    
        end
        
        
    
        
        pyramids = getMaxPyramids(currImage, baseImage);
    
        %Create vector for iterations, reducing iterations at
        %higher levels (which are the greatest time burden)
        %improves runtime without much effect on performance
        iterations = round(linspace(numIterations, 10, pyramids));
    
        %Perform non-rigid registration using the previous
        %registered image, allows transformations to combine so all
        %images align
        
        [D,currReg] = imregdemons(currImage, baseImage, iterations,...
    'PyramidLevels', pyramids, 'AccumulatedFieldSmoothing', smoothing, 'DisplayWaitbar', false);
    
        displacementFields{selected} = D;
    
    
        registeredBulkI{selected} = currReg;

        
    
        %Adjust centers so all align
        centersXRounded = round(centers{selected}(:,1));
        centersYRounded = round(centers{selected}(:,2));
    
    
    
    
        adjustedX = zeros(1, length(centersXRounded));
        adjustedY = zeros(1, length(centersYRounded));
        
        sourceX(:,:) = meshX(:,:) + D(:,:,1);
        sourceY(:,:) = meshY(:,:) + D(:,:,2);
         

        
        
        %Apply the displacement field to the centers to align
        %them for comparison
        for center = 1:length(centersXRounded)
            diffFromSource = abs(sourceX-centersXRounded(center)).^2 + abs(sourceY-centersYRounded(center)).^2;
            [~,closestIndex] = min(diffFromSource(:),[],'all',"linear");
            [y, x] = ind2sub(size(D), closestIndex);
            adjustedX(center) = x;
            adjustedY(center) = y;
  
            
        end
    
    
        
    
        
        
    
        adjustedCenters{selected} = [adjustedX' adjustedY'];
    end
end


function [adjustedCenters, totalTform] = rigidRegister(cc, centers, bulkI, regMethod, rotate)
    adjustedCenters = cell(1,length(bulkI));
    totalTform = cell(1,length(bulkI));
    registeredBulkI = cell(1,length(bulkI));
    adjustedCenters{1} = centers{1};
    firstImageCCLabel = labelmatrix(cc{1});
    if strcmp(regMethod, "Mask")
        registeredBulkI{1} = logical(firstImageCCLabel);
        
        
    elseif strcmp(regMethod, "Image")
   
        registeredBulkI{1} = bulkI{1};
        %app.RegisteredBulkI{app.1} = imbinarize(app.BulkI{app.1}, 'adaptive', 'ForegroundPolarity', 'Dark', 'Sensitivity', 0.6);
    
        
    else
        tempBulkI = zeros(size(firstImageCCLabel));
        for center = 1:size(centers{1}, 1)
            tempBulkI(round(centers{1}(center, 2)), round(centers{1}(center, 1))) = 1;
        end
        registeredBulkI{1} = imdilate(tempBulkI, ones(5));
       

    end
    totalTform{1} = affine2d(eye(3));
    firstPerim = bwperim(firstImageCCLabel);
    firstPerim = imdilate(firstPerim, ones(3,3));
    registeredPerimeter{1} = firstPerim;
    
   
    
    
    for selected = 2:length(bulkI)
        prev = selected - 1;
    %                 
    
        %firstImageCCLabel = labelmatrix(app.CC{prev});
        secondImageCCLabel = labelmatrix(cc{selected});
        
        %Perform registration with the previously registered image
        %so that displacement fields for all images align
        
        baseImage = registeredBulkI{prev};
        
        if strcmp(regMethod, "Mask")
            currImage = logical(secondImageCCLabel);
            
            
        elseif strcmp(regMethod, "Image")
       
            currImage = bulkI{selected};
            %app.RegisteredBulkI{app.1} = imbinarize(app.BulkI{app.1}, 'adaptive', 'ForegroundPolarity', 'Dark', 'Sensitivity', 0.6);
        
            
        else
            tempBulkI = zeros(size(secondImageCCLabel));
            for center = 1:size(centers{selected}, 1)
                tempBulkI(round(centers{selected}(center, 2)), round(centers{selected}(center, 1))) = 1;
            end
            currImage = imdilate(tempBulkI, ones(3));
            
    
        end
        
      
        
    
        if ~rotate
            [totalTformCurr, peakcorr] = imregcorr(currImage, baseImage, 'translation');
        else
            [totalTformCurr, peakcorr] = imregcorr(currImage, baseImage, 'rigid');
        end
    
     
        
    
        totalTform{selected} = totalTformCurr;
        sameAsInput = affineOutputView(size(baseImage), totalTformCurr, "BoundsStyle","SameAsInput");
        registeredBulkI{selected} = imwarp(currImage, totalTform{selected},"OutputView",sameAsInput);
              
        
        
    

    
        %Adjust centers so all align
        centersXRounded = round(centers{selected}(:,1));
        centersYRounded = round(centers{selected}(:,2));
    
    
    
    
        adjustedX = zeros(1, length(centersXRounded));
        adjustedY = zeros(1, length(centersYRounded));
        
       
         
   
        
        
        %Apply the displacement field to the centers to align
        %them for comparison
        for center = 1:length(centersXRounded)
            oldCoord = [centersXRounded(center) centersYRounded(center) 1];
            newCoord = oldCoord*totalTformCurr.T;
            
            adjustedX(center) = newCoord(1);
            adjustedY(center) = newCoord(2);
   
            
        end
    
    
        
    
        
        
        adjustedCenters{selected} = [adjustedX' adjustedY'];
    end
end


function maxPyr = getMaxPyramids(curr, base)
    dims = [size(curr), size(base)];
    minDim = min(dims);
    i = 1;
    
    while 2^i <= minDim
        i = i + 1;                
    end
    maxPyr = i - 1;
end


function matchingIDs = trackOrganoids(adjustedCenters, maxDistance)
    numObjects = zeros(1, length(adjustedCenters));
    for selected = 1:length(adjustedCenters)
        numObjects(selected) = size(adjustedCenters{selected}, 1);
    end
    
    matchingIDs = cell(1, length(adjustedCenters));
    
    
    [xCoordsAdj, yCoordsAdj] = getMetrics(adjustedCenters, numObjects);
    numNotAssignedPerImage = cell(1, length(adjustedCenters) - 1);
    
   
        matchingIDs{1} = 1:size(adjustedCenters{1},1);
       
        
        
        for selected = 2:length(adjustedCenters)
            
            
            prev = selected - 1;

            while prev >= 1
                
                xs = [xCoordsAdj{prev}' xCoordsAdj{selected}'];
                ys = [yCoordsAdj{prev}' yCoordsAdj{selected}'];

                xDist = findDifference(xs, selected, prev, numObjects);
                
                yDist = findDifference(ys, selected, prev, numObjects);
                xyDist = sqrt(xDist.^2+yDist.^2);
       

                

                distMatrix = xyDist;
                
               

                distMatrix(xyDist > (maxDistance)) = nan;
                
                

                
                prevMatchingID = matchingIDs{prev};

                %If this is the first iteration backwards, allocate
                %otherwise, clear out the values that have already been
                %matched
                if prev == selected - 1
                    
                   

                    matchingID = zeros(1, size(adjustedCenters{selected},1));
                else    
                    %Clear out spheroids from this image that have been
                    %matched
                    distMatrix(:, matchingID ~= 0) = nan;
                    alreadyUsedID = ismember(prevMatchingID, matchingID);
                    %Clear out spheroids from previous image that were
                    %matched to this image already
                    distMatrix(alreadyUsedID, :) = nan;
                end
                
                %Check if all spheroids have been assigned
                allAssigned = false;
                if nnz(~isnan(distMatrix)) == 0
                        allAssigned = true;
                end

                %Perform assignment for the current iteration of image
                %and a previous image
                while ~allAssigned

                    temp = distMatrix;
                    %temp(isnan(distMatrix())) = maxVal;
                    
                    %Find index of minimum value in distance matrix
                    [~,X] = min(temp(:));
                    [R,C] = ind2sub(size(distMatrix),X);
                    
                    
                    %Assign to matching ID, and then clear in current
                    %and prev so that they cannot be assigned again
                    matchingID(C) = prevMatchingID(R);

                    distMatrix(R,:) = nan;
                    distMatrix(:,C) = nan;
                    
                    
                                       

                    if nnz(~isnan(distMatrix)) == 0
                        allAssigned = true;
                    end
                    
                end
                prev = prev - 1;
            end

            
            %Give IDs to all spheroids from this current image that were not
            %assigned 
            numNotAssigned = 0;
            nextVal = getMaxMatchingID(matchingIDs) + 1;
            for i = 1:length(matchingID)
                if matchingID(i) == 0
                    numNotAssigned = numNotAssigned + 1;
                    matchingID(i) = nextVal;
                    nextVal = nextVal + 1;
                end
            end
            numNotAssignedPerImage{selected-1} = [numNotAssigned, size(adjustedCenters{selected},1)];


            matchingIDs{selected} = matchingID;
            
            
            
        end

       

        
   
end


function [xCoordsAdj, yCoordsAdj] = getMetrics(adjustedCenters, numObjects)
            
            xCoordsAdj = cell(1, length(adjustedCenters));
            yCoordsAdj = cell(1, length(adjustedCenters));
            
          
            maxNumObjects = max(numObjects);

            
            for selected = 1:length(adjustedCenters)
                numObjectsCurr = numObjects(selected);
              
                xcoordsAdjCurr = zeros(1, maxNumObjects);
                ycoordsAdjCurr = zeros(1, maxNumObjects);
                

               
                

                

               




                

                
                
             
                xcoordsAdjCurr(1:numObjectsCurr) = adjustedCenters{selected}(:,1);
                ycoordsAdjCurr(1:numObjectsCurr) = adjustedCenters{selected}(:,2);
                xCoordsAdj{selected} = xcoordsAdjCurr;
                yCoordsAdj{selected} = ycoordsAdjCurr;

                
                

                





                
            end
            
            
end


function maxNum = getMaxMatchingID(matchingIDs)
    maxNum = 0;
    
    for i = 1:length(matchingIDs)
        tempMax = max(matchingIDs{i});
        if tempMax > maxNum
            maxNum = tempMax;
        end
    end
    
end


function [correctComparisons, totalComparisons, correctRetained, totalRetained] = evaluatePerformance(matchingIDs, manualTrack)
    totalComparisons = zeros(1,length(matchingIDs));
    correctComparisons = zeros(1,length(matchingIDs));
    totalRetained = zeros(1,length(matchingIDs));
    correctRetained = zeros(1,length(matchingIDs));
    totalRetained(1) = sum(~isnan(manualTrack(:,1)));
    correctRetained(1) = sum(~isnan(manualTrack(:,1)));
    correspondingIDs = getCorrespondingIDs(matchingIDs, manualTrack);
    for image = 2:length(matchingIDs) 
        totalRetainedCurr = 0;
        correctRetainedCurr = 0;
        totalComparisonsCurr = 0;
        correctComparisonsCurr = 0;
        for track = 1:size(manualTrack, 1)
            track
            if ~isnan(manualTrack(track,image-1))
                totalComparisonsCurr = totalComparisonsCurr + 1;
                'comparison'
                if isnan(manualTrack(track,image))
                    
                    if isempty(find(matchingIDs{image} == correspondingIDs(track)))
                        
                        correctComparisonsCurr = correctComparisonsCurr + 1;
                        'correct empty'
                    else
                        'incorrect'
                    end

                else
                  
                    
                    

                    if matchingIDs{image-1}(manualTrack(track, image-1)) == matchingIDs{image}(manualTrack(track, image))
                        correctComparisonsCurr = correctComparisonsCurr + 1;
                    else
                        'incorrect num'
                    end

                end

            
            end
            
            if ~isnan(manualTrack(track,image))
                totalRetainedCurr = totalRetainedCurr + 1;
                
                
                if matchingIDs{image}(manualTrack(track, image)) == correspondingIDs(track)
                    correctRetainedCurr = correctRetainedCurr + 1;

                end
            end
            
        end
        totalRetained(image) = totalRetainedCurr;
        correctRetained(image) = correctRetainedCurr;
        totalComparisons(image) = totalComparisonsCurr;
        correctComparisons(image) = correctComparisonsCurr;
    end
end

function correspondingIDs = getCorrespondingIDs(matchingIDs, manualTrack)
  
    correspondingIDs = zeros(1, size(manualTrack, 1));
    for track = 1:size(manualTrack, 1)
        firstNonNan = find(~isnan(manualTrack(track,:)), 1, "first");
        firstNonNan = firstNonNan(1);
        correspondingID = (matchingIDs{firstNonNan}(manualTrack(track, firstNonNan)));
        correspondingIDs(track) = correspondingID(1);
    end
end


function differenceMatrix = findDifference(observations, current, prev, numObjects)
    
    numObservationsFirst = numObjects(prev);
    numObservationsSecond = numObjects(current);
    differenceMatrix = zeros(numObservationsFirst, numObservationsSecond);
    for i = 1:numObservationsFirst
        for j = 1:numObservationsSecond
            differenceMatrix(i,j) = abs(observations(i,1) - observations(j,2));
        end
    end
end



