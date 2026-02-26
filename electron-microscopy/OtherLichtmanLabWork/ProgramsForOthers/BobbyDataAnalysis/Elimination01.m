
while 1
clear all



jSiz = [ 7 7 ];
jNum = 9; 
jNum = fix(sqrt(jNum))^2;
aNum = 5;
colNum = uint8(aNum + 1);
colormap lines(100)

extRate = .5;
resource = 10;
occReward = .0001;
sizPunish = .4;

reps = 1000;
showEvery = 10;

%% make junctions
juncs = fix(rand(jSiz(1),jSiz(2),jNum)*aNum)+1;
buffJuncs = zeros(jSiz(1)+2,jSiz(2)+2,jNum);

%% make axons
for a = 1:aNum
   eRat(a) = extRate;
   nativeRetract(a) = .1;
   
end

%% iterate elemination
clear occ
for r = 1: reps
    
    for j = 1: jNum  %get occupancy
        junc = juncs(:,:,j);
        occ(j,:) = hist(juncs(:),1:aNum);
    end
    occ = occ/prod(jSiz);  %get junction occupancy
    aSiz = hist(juncs(:),1:aNum);  %get arbor sizes
    aSiz = aSiz/sum(aSiz);
    
    axRet = nativeRetract + (sizPunish .* aSiz); % axon specific retraction  
    
    for j = 1 : jNum
        probRet = axRet - (occ(j,:) .* occReward); 
        randRet = rand(jSiz);
        juncRet = zeros(jSiz);
        
        for a = 1:aNum  %  make retraction matrix
           juncRet(juncs(:,:,j)==a)=probRet(a);            
        end
        junc = juncs(:,:,j);    
        junc(juncRet>randRet) = 0; %RETRACT
        juncs(:,:,j) = junc;
        
    end
    
    %% Run extensions
    empties = find(~juncs);
    [eY eX eZ] = ind2sub([jSiz(1) jSiz(2) jNum], empties);
    eYX = [eY eX];
    shift = [-1 0; 0 -1; 1 0; 0 1];
    
    buffJuncs(2:size(juncs,1)+1,2:size(juncs,2)+1,:)=juncs;
    for e = 1:length(empties)
        sur(:,1) = shift(:,1) + eY(e)+1;
        sur(:,2) = shift(:,2) + eX(e)+1;
        sur(:,3) = eZ(e);
        targA = sub2ind([jSiz(1)+2 jSiz(2)+2 jNum],sur(:,1),sur(:,2),sur(:,3)); 
        surA = buffJuncs(targA);
        nearby = find(surA);
        if ~isempty(nearby) %if there is a nearby axon
            useAx = surA(nearby(fix(rand*length(nearby)+1))); %pick surround axon
            juncs(eY(e),eX(e),eZ(e)) = useAx;
        end
        
        
    end
    
    
    %% Display
    if mod(r,showEvery) ==0
        subNum = fix(sqrt(jNum));
        for j = 1: jNum
            subplot(subNum,subNum,j)
            image(juncs(:,:,j))
        end
        pause(.01)

        r
    end
end


end

