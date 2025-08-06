

if exist('mandersSegInGoldAvgPerim', 'var') == 0
    compareToGoldCC;
end


numSegments  = 0;
for image = 1:length(mandersSegInGoldTotalPerim)
    numSegments = numSegments + length(mandersSegInGoldTotalPerim{image});
end

mandersSegInGoldPerimTogether = zeros(4, numSegments, numDilationLevels);



for condition = 1:4
    startIndex = 0;
    for image = 1:length(mandersSegInGoldTotalPerim)
        
        
        mandersSegInGoldPerimTogether(condition, startIndex + 1: startIndex+length(mandersSegInGoldTotalPerim{image}), :) = mandersSegInGoldTotalPerim{image}(condition, :, :);
        
        
        startIndex = startIndex + length(mandersSegInGoldTotalPerim{image});
    end
end
rows = 1;
















figure(2)
clf


i = 1;
    sortedCdf1 = sort(mandersSegInGoldPerimTogether(1,:,i));
    sortedCdf2 = sort(mandersSegInGoldPerimTogether(2,:,i));
    mandersValues1 = [.7368 .8085 .8052 .5709 .7116 .5878];
    percentiles1 = zeros(1,length(mandersValues1));
    for j = 1:length(percentiles1)
       

        percentiles1(j) = findPercentile(sortedCdf1, mandersValues1(j));
    end
    mandersValues2 = [.8803 .8598 .9702 .7594 .9326 1];
    percentiles2 = zeros(1,length(mandersValues2));
    for j = 1:length(percentiles2)
        percentiles2(j) = findPercentile(sortedCdf2, mandersValues2(j));
    end

    h = cdfplot(mandersSegInGoldPerimTogether(1,:,i));
    h.Color = "black";
    h.LineStyle = '--';
    h.LineWidth = 1;
    
   



  
    colors = {[55 126 184]/255, [77 175 74]/255, [255 192 203]/255, [255 127 0]/255, [166 86 40]/255, [152 78 163]/255};

    hold on
    h = cdfplot(mandersSegInGoldPerimTogether(2,:,i));
    h.Color = "black";
    
    
    hold on
    for point = 1:length(mandersValues1)
        
        plot(mandersValues1(point), percentiles1(point), '.', 'Color', colors{point}, 'MarkerSize', 20);
        plot(mandersValues2(point), percentiles2(point), '.', 'Color', colors{point}, 'MarkerSize', 20);
    end
    
    
    hold on

    
    xlim([0.5,1])
    ylim([0 0.7])
    hold on
    legend('OrganoSeg', 'OrganoSeg2', 'Location', 'northwest')
    dilationString = strcat('Dilation=', num2str(dilations(i)));
    title("CDF of Mander's Coefficient")
    
    xlabel("Mander's Colocalization Coefficient")
    hold off
    
    [h,p] = kstest2(mandersSegInGoldPerimTogether(1,:,1), mandersSegInGoldPerimTogether(2,:,1))
    subtitle(sprintf("P = %.5g", p), 'FontSize', 9)
    pbaspect([1,1,1])

    


hold off






%% 



function percentile = findPercentile(sortedcdf, mander)
index = find(round(sortedcdf,4) == mander, 1);
percentile = index/length(sortedcdf(~isnan(sortedcdf)));
end