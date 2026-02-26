function[kern] = jmkern(type,var1,var2,var3)

if ~exist('var1','var')
    var1 = 10;
end

if strcmp(type,'ball')
    
    ballRad = round(var1/2);
    
    ball = ones(ceil(ballRad)*2+1,ceil(ballRad)*2+1,ceil(ballRad)*2+1);
    [y x z] = ind2sub(size(ball),find(ball));
    dists = sqrt((y-mean(y)).^2+(x-mean(x)).^2+ (z-mean(z)).^2);
    ball(:) = dists;
    
    ball2 = ball<=ballRad;
    ballSum = sum(ball2,3);
    ballSum = ballSum/max(ballSum(:));
    ballSum(ballSum<.5) = 0;
    kern = ballSum;
elseif strcmp(type, 'disk')
    
      ballRad = round(var1/2);
    
    ball = ones(ceil(ballRad)*2+1,ceil(ballRad)*2+1,ceil(ballRad)*2+1);
    [y x z] = ind2sub(size(ball),find(ball));
    dists = sqrt((y-mean(y)).^2+(x-mean(x)).^2+ (z-mean(z)).^2);
    ball(:) = dists;
    
    ball2 = ball<=ballRad;
    ballSum = sum(ball2,3);
    ballSum = ballSum/max(ballSum(:));
    ballSum(ballSum<.5) = 0;
    kern = ballSum;
    kern = kern>0;
    
    
elseif strcmp(type, 'ring')
    
    ballRad = round(var1/2);
    if ~exist('var2','var')
        var2 = ceil(var1/3);
    end
    ringWidth = var2;
    
    ball = ones(ceil(ballRad)*2+1,ceil(ballRad)*2+1);
    [y x] = ind2sub(size(ball),find(ball));
    dists = sqrt((y-mean(y)).^2+(x-mean(x)).^2);
    ball(:) = dists;
    
    ball2 = ball<=(ballRad);
    ball2(ball<(ballRad-ringWidth)) = 0;
    ballSum = ball2;
    kern = ballSum;
    
elseif strcmp(type, 'bar')
    
    
    if ~exist('var2','var')
        var2 = ceil(var1/3);
    end
    kern = ones(var1,var2);
    
    
elseif strcmp(type, 'I')
    
    if ~exist('var2', 'var')
        var2 = round(var1/3);
    end
    
    kern = [
        1 1 1 1 ;....
        0 1 1 0 ;....
        0 1 1 0 ;....
        0 1 1 0 ;....
        1 1 1 1];
    
    kern = imresize(kern,[var1 var2],'nearest');
    
elseif strcmp(lower(type),'x')
    
   if ~exist('var2', 'var')
        var2 = round(var1/3);
   end
    
    
   
   kern = [ 1 1 0 0 0 0 0 0 0 1 1 ; ...
       1 1 1 0 0 0 0 0 1 1 1 ; ...
       0 1 1 1 0 0 0 1 1 1 0 ; ...
       0 0 1 1 1 0 1 1 1 0 0 ; ...
       0 0 0 1 1 1 1 1 0 0 0 ; ...
       0 0 0 0 1 1 1 0 0 0 0 ; ...
       0 0 0 1 1 1 1 1 0 0 0 ; ...
       0 0 1 1 1 0 1 1 1 0 0 ; ...
       0 1 1 1 0 0 0 1 1 1 0 ; ...
       1 1 1 0 0 0 0 0 1 1 1 ; ...
       1 1 0 0 0 0 0 0 0 1 1 ];
   
   kern = imresize(kern,[var1 var2],'nearest');
   
   
end





kern = double(kern)/max(kern(:));

















