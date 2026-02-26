%%Extract contrast from one xml file and apply it to all


%% Get xml file to use and find all others.
fclose('all')
[TFN, TPN] = uigetfile('*.*');  % Get file and directory


%% Get contrast from targeted xml

fid = fopen([TPN TFN])
file = textscan(fid, '%s', 'delimiter', '\n','whitespace', '');
lines = file{1};
fclose(fid)

%%Get Contrast and Brightness
for i = 1:length(lines)
    line = lines{i};
    conPos = findstr(line,'contrast'); % get positions of contrast
    if ~isempty(conPos)
        quotePos = findstr(line,'"');  % get positions of quotes
        getPos = quotePos(quotePos>conPos); % look at qotes after 'contrast'
        useContrast = line(getPos(1)+1:getPos(2)-1) % grab contrast
        useBrightness = line(getPos(3)+1:getPos(4)-1)  % grab brightness
        break
    end
end

%pause 

%% find xml to change
dots = find(TFN == '.');
GFN = TFN(1:dots(length(dots)));  %extract generic file name
namLength = length(GFN);

xList = {};  %%initialize list of xml files

TPNd = dir(TPN);
for i = 1: length(TPNd)

    nam = TPNd(i).name;
    if length(nam) == length(TFN);
        if strcmp(nam(1:namLength), GFN);
            xList{length(xList)+1} = nam;
        end % end if the right name
    end  % end if right length
end  % end looking for xml files


%% Enter new brightness and contrast

for x = 1: length(xList)
    percentDone = x/length(xList) * 100
    
    xFN = xList{x};
    fid = fopen([TPN xFN]);
    file = textscan(fid, '%s', 'delimiter', '\n','whitespace', '');
    fclose(fid);
    lines = file{1};

    %%Insert new contrast and brightness
    for i = 1:length(lines)
        line = lines{i};
        conPos = findstr(line,'contrast'); % get positions of contrast
        if ~isempty(conPos)
            quotePos = findstr(line,'"');  % get positions of quotes
            getPos = quotePos(quotePos>conPos); % look at qotes after 'contrast'
            contrast = line(getPos(1):getPos(2)); % grab contrast
            rightness = line(getPos(3):getPos(4));  % grab brightness
            newLine = [line(1:getPos(1)) useContrast line(getPos(2):getPos(3)) useBrightness line(getPos(4):length(line))];
            lines{i} = newLine;
            break
        end
    end

    %%Convert to matrix for writing
    
    maxlength = 145;  %% get number of columns
    for i = 1:length(lines)
        maxlength = max(maxlength,length(lines{i}));
    end
    myText = char(zeros(length(lines),maxlength)+32);  %initialize text file

    for i = 1:length(lines) %write lines in to text variable
        lineText = lines{i};
        if length(lineText)
            myText(i,1:length(lineText)) = lineText;
        end
    end

    %%Write text into xml file
    dlmwrite([TPN xFN],myText,'delimiter','','newline','pc');

end % run all xml files



