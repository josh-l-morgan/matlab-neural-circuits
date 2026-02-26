function[badTargNum] = testReciprocal(con2);

%%

%con2 = vox.conMat;
sourceCon = con2*0;
[checkY checkX] = find(con2>0);

sourceRef = zeros(1,26);

%sourceRef =     [6     5     4     3     2     1    18    17    16    15    14    13    12    11    10     9     8     7    26    25    24    23    22    21    20    19];

for i = 1:length(checkY)
    val = con2(checkY(i),checkX(i));
    
    %get source position
    if sourceRef(checkX(i))==0
        sourcePos = find(con2(val,:)==checkY(i))  % find source Position
        if ~isempty(sourcePos) % if good value
            sourceRef(checkX(i)) = sourcePos
        end
    else
        sourcePos = sourceRef(checkX(i)); %retrieve source position
    end
    
    % record error
    if isempty(sourcePos) %no position found
        sourceCon(checkY(i),checkX(i)) = 0; % if couldnt find recipricol
    else
        sourceVal = con2(val,sourcePos); 
        if sourceVal ~= checkY(i) %position contains incorrect number
            sourceCon(checkY(i),checkX(i)) = sourceVal; % if couldnt find recipricol
        end
    end
    
end



mistakeNum = sum(sourceCon(:)>0);
mistakeID = unique(sourceCon(sourceCon>0));
badTargNum = length(mistakeID);


