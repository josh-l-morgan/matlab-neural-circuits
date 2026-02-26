function[sqrC] = squaredCluster(testMat)

testMat = testMat(:);
%sqrC = sqrt(mean(testMat(:).*testMat(:)));

%% sqrt(sum(testMat))/sum(testMat)    : 1
%sqrC = sqrt(sum(testMat .* testMat))/sum(testMat);



%% 1/sum(testMat)    : 1
% sqrC = sum(testMat.*testMat)/sum(abs(testMat))^2;



%% 1    : max cluster size %% JM synaptically weighted Average of cluster size
%sqrC = sum(testMat.*testMat)/sum(abs(testMat));


%% Victor normalized: 
% sqrC = (sum(testMat))^2 / (sum(testMat^2)* lentgh(testMat))  

%% Rms error ratio
     
%% just RMS error

%sqrC = rms(testMat-meanCluster)/meanCluster;


%% I dont know
%sqrC = sum(testMat.^2)/sum(abs(testMat))/meanCluster;

%% Qartile

% sortMat = sort(testMat,'descend');
% qrt = ceil(length(sortMat)/10);
% 
% sqrC = (mean(sortMat(1:qrt))-mean(sortMat(end-qrt+1:end)))/mean(sortMat);
% sqrC = (mean(sortMat(1:qrt))-mean(sortMat))/mean(sortMat);

%% rms normal

%sqrC = rms(testMat-mean(testMat))/mean(testMat);


%% skewness
%sqrC = skewness(testMat);


%% normalized most common cluster size
%sqrC = (sum(testMat.^2)/sum(testMat))/mean(testMat);


%% 

%sqrC = sqrt(sum(testMat.^2))/mean(testMat);

%%

sqrC = rms((testMat-mean(testMat))/mean(testMat));


%%

%sqrC = (mean(((testMat-mean(testMat))/mean(testMat)).^4))^(1/4);




