function[sqrC] = squaredCuster(testMat)

testMat = testMat(:);
%sqrC = sqrt(mean(testMat(:).*testMat(:)));

%% sqrt(sum(testMat))/sum(testMat)    : 1
%sqrC = sqrt(sum(testMat .* testMat))/sum(testMat);



%% 1/sum(testMat)    : 1
sqrC = sum(testMat.*testMat)/sum(testMat)^2;


