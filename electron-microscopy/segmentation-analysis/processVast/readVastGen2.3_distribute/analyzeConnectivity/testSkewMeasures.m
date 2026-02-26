

%% Test Cluster index

testDists = {[2 2 2 2] [3 2 1 2] [6 2 0 0] [6 2 0 0 0 0 ] [8 0 0 0] [8 0 0 0 0 0 0 0 ]}

clear cIs 

subplot(1,1,1)

for i = 1 :length(testDists)
    
       tD = testDists{i};

    clear clusterIndex
    meanCluster = sum(tD)/length(tD);
difCluster = tD-meanCluster;
mostCluster = sum(tD) - meanCluster + meanCluster * (length(tD)-1);
mostCluster = sum(tD)  + meanCluster * (length(tD)-2);

%mostClusterSqr = (sum(tD)-meanCluster)^2 + (meanCluster^2) * (length(tD)-1);

leastCluster = 0;

skewMat = tD;
skewMat(1) = sum(tD);
mostClusterSqr = (sum((skewMat-meanCluster).^2));
mostClusterSqr2 = sqrt(mean((skewMat-meanCluster).^2));

    
   %clusterIndex(1) = sum(tD.^2)/sum(tD); %good - weighted average
   %clusterIndex(2) = sum(tD.^2)/sum(tD)^2; % good - cluster index
   %clusterIndex(3) = sqrt(sum(tD.^2))/sum(tD); % weak compared to #2
   %clusterIndex(4) = (sum(tD))^2 / (sum(tD.^2)* length(tD)) 
   
  meanTD = mean(tD);
  medianTD = median(tD);
  stdTD = std(tD);
  n = length(tD);
  
  clusterIndex(1) =n/((n-1)*(n-2)) * sum(((tD - mean(tD))/std(tD)).^3); %adjusted Fisher-Pearson standardized moment coefficient
  col{1} = 'r';
   clusterIndex(2) =   (mean(tD) - median(tD))/std(tD); %non parametric skew
   col{2} = 'g';
   clusterIndex(3) = rms(tD-mean(tD))/mean(tD); %same as above
   col{3} = 'b';
   clusterIndex(4) =   rms((tD-mean(tD))/mean(tD));
   col{4} = 'c'; 
   clusterIndex(5) = skewness(tD,0)+.02;
   col{5} = 'y'; 
   clusterIndex(6) = (sum(tD.^2)/sum(tD))/mean(tD);
   col{6} = 'm';
   clusterIndex(7) = sqrt(mean((tD - mean(tD)).^2))/mean(tD)+.03;
   col{7} = 'k';
   
   cIs(:,i) = clusterIndex;
    
end



for i = 1:size(cIs,1)
    
    try useCol = col{i};
    catch err
        useCol = 'k';
    end
    plot(cIs(i,:),useCol)
    hold on
end
hold off
    






