

%% Test Cluster index

testDists = {[2 2 2 2] [3 2 1 2] [6 2 0 0] [6 2 0 0 0 0 ] [8 0 0 0] [8 0 0 0 0 0 0 0 ]}

clear cIs 



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
  clusterIndex(1) = sum((tD-meanCluster).^2)/mostClusterSqr;
   clusterIndex(2) =   sqrt(sum((tD-meanCluster).^2))/sqrt(mostClusterSqr);
   clusterIndex(3) = rms(tD-meanCluster)/rms(skewMat-meanCluster); %same as above
  % clusterIndex(4) =   sqrt(mean((tD-meanCluster).^2))/sqrt(mostClusterSqr2);
   clusterIndex(4) =   sqrt(mean((tD-meanCluster).^2))/sqrt(mostClusterSqr2);
    clusterIndex(5) = rms(tD-meanCluster)/meanCluster;

   cIs(:,i) = clusterIndex;
    
end


plot(cIs')
%%


(td - sum(tD)/length(tD))^2/sum(tD)^2


sum(tD)^2-(sum(tD)/length(tD))




