function[ballSum] = ballPic(ballRad)


ball = ones(ceil(ballRad)*2+1,ceil(ballRad)*2+1,ceil(ballRad)*2+1);
[y x z] = ind2sub(size(ball),find(ball));
dists = sqrt((y-mean(y)).^2+(x-mean(x)).^2+ (z-mean(z)).^2);
ball(:) = dists;

ball2 = ball<=ballRad;
ballSum = sum(ball2,3);
ballSum = ballSum/max(ballSum(:));
ballSum(ballSum<.5) = 0;

% 
% ball1 = ball;
% ball1(ball1>=ballRad) = 0;
% ballSum = max(ball1>0,[],3);
% % 
% % ball2 = ball;
% % ball2(ball2>5) = 0;
% ballSum2 = sum(ball2>0,3);
% ballSum2 = ballSum2 * 200/max(ballSum2(:));
% 
% ballSum = ballSum1/100+ballSum2;