function[bSub] = ballz(anchors, ballRad);
%%function[bSub] = ballz(anchors, balRad);


        
        
        %%Ball
        if ~exist('ballRad','var')
            ballRad = 6;
        end
        
        ball = ones(ballRad*2+1,ballRad*2+1,ballRad*2+1);
        [y x z] = ind2sub(size(ball),find(ball));
        dists = sqrt((y-mean(y)).^2+(x-mean(x)).^2+ (z-mean(z)).^2);
        ball(:) = dists;
%         ball(ball>ballRad) = 0;
%         ball = ball>0;
        ballInd = find(ball<=ballRad);
        [y x z] = ind2sub(size(ball),ballInd);
        yxz = [y x z]- mean(y);
        
        bSize = size(yxz,1);
        %bSub = zeros(size(anchors,1) *bSize,3);
        
        for i = 1:size(anchors,1);
            %bSub((i-1)* bSize+1 : (i-1)*bSize+bSize,1:3) = ...
            %    repmat(anchors(i,:),[bSize 1]) + yxz;
            bSub{i} = repmat(anchors(i,:),[bSize 1]) + yxz;
        end
        
        
        