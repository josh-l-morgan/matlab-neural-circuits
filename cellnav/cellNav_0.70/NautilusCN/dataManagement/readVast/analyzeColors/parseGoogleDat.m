function[dat] = parseGoogleDat(dat)


if ~exist('dat','var')    
    load('MPN.mat')
    load([MPN 'dat.mat'])
end

googleLinks = dat.googleLinks;
clear dat
dat.googleLinks = googleLinks;

rawDat = GetGoogleSpreadsheet2(googleLinks.datSheet,googleLinks.datGID);
rawAlias = GetGoogleSpreadsheet2(googleLinks.datSheet,googleLinks.aliasGID);


%% parse alias
c = 0;
alias = {};
for y = 1:size(rawAlias,1)
    
    ids = [];
    for x = 1:size(rawAlias,2)
        str = rawAlias{y,x};        
        if ~isempty(str)
            num = sscanf(str,'%d');
            ids = [ids num];
        end
    end
    
    if ~isempty(ids)
        c = c+1;
        alias{c} = ids;
    end
end

dat.alias = alias;




%% parse cellDat

dNum = size(rawDat,1)-1;

for x = 1:size(rawDat,2);
   
        str = rawDat{1,x};
        if isempty(str)
            str = sprintf('col %d',x);
        end
        dat.columnNames{x} = str;
end

dat.rawDat = rawDat;

%% Parse IDs
dat.cid = zeros(dNum,1);
for y = 2:size(rawDat,1)
   
    str = rawDat{y,1};
    nums = sscanf(str,'%d');
    if ~isempty(nums)
        dat.cid(y-1) = nums(1);
    end

end


%% parse cell type
dat.typeNames = {'rgc' 'tcr' 'lin' 'unk','phot','hrz','bpc','amc','mul'};
dat.type = zeros(dNum,1);
for y = 2:size(rawDat,1)
    
   dat.type(y-1) = 0; 
   
   str = lower(rawDat{y,2});
   notSpc = find(~(str == ' '));
   if isempty(str), str = ' ';
   else
    str = str(notSpc(1) : notSpc(end));
   end
       
   
   targ = find(strcmp(dat.typeNames,str),1);
   if ~isempty(targ)
    dat.type(y-1) = targ;
   else
       dat.typeNames = cat(2,dat.typeNames,{str});
       dat.type(y-1) = length(dat.typeNames);
   end
           
end


%% parse on/off
dat.onOffNames = {'on' 'off'};
dat.onOff = zeros(dNum,length(dat.onOffNames));
for y = 2:size(rawDat,1)
    
   dat.onOff(y-1,1:2) = 0; 
   str = lower(rawDat{y,3});
   for t = 1:length(dat.onOffNames)
       %str,dat.onOffNames{t},pause
       if sum(regexp(str,dat.onOffNames{t}))
         dat.onOff(y-1,t) = 1;
       end
   end          
end



%% parse subType
for i = 1:length(dat.typeNames),dat.subTypeNames{i} = {};end
dat.subTypeNames{8} = {'vgc'};
dat.subTypeNames{7} = {'bc1' 'bc2' 'bc3a' 'bc3b' 'bc4' 'bc5i' 'bc5o' 'bc5t' 'bc6' 'bc7' 'bc8' 'bc11' 'unk','on','off'};
dat.subTypeNames{1} = {'unk' 'on' 'off' 'on/off' '1ni' '1no' '1wt' '25' '27' '28'  '2an' '2aw' '2i' '2o' ...
    '37c' '37d' '37r' '37v' '3i' '3o' '4i' '4on' '4ow' '51' '5si' '5so' '5ti' ...
    '5to' '63' '5sn' '5sw' '6t' '72' '73' '7id' '7ir' '7iv' '7o' '81i' '81o' '82n' ...
    '82wi' '82wo' '85' '8n' '8w' '91' '915' '9n' '9w' '37' '4' '5' '6' '7'};
dat.subType = zeros(dNum,1);
for y = 2:size(rawDat,1)
    
    str = lower(rawDat{y,4});
     notSpc = find(~(str == ' '));
   if isempty(str), str = ' ';
   else
    str = str(notSpc(1) : notSpc(end));
   end
   
    if isempty(str)
        dat.subType(y-1) = 0;
    elseif dat.type(y-1) == 0;
        dat.subType(y-1) = 0;
    else
        typeID = dat.type(y-1);
        checkNames = dat.subTypeNames{typeID};
        
        isName = strcmp(checkNames,str);
        if sum(isName)
            dat.subType(y-1) = find(isName,1);
        else
           
            %%Make new name
            dat.subTypeNames{typeID} = cat(2,dat.subTypeNames{typeID},str);
            dat.subType(y-1) = length(checkNames) + 1;
        end
            
        
    end
end




%% parse tracing
dat.tracingNames = {'none' 'fragment' 'partial' 'full' 'error'};
dat.tracing = zeros(dNum,1);
for y = 2:size(rawDat,1)
   dat.tracing(y-1) = 0; 
   str = lower(rawDat{y,5});
   for t = 1:length(dat.tracingNames)
       if sum(regexp(str,dat.tracingNames{t}))
         dat.tracing(y-1) = t;
         break
       end
   end          
end




%% parse seed Position
dat.seedPosition = zeros(dNum,3);
for y = 2:size(rawDat,1)
    
    dat.seedPosition(y-1,1:3) = 0;
    str = rawDat{y,6};
    for t = 1:length(dat.subTypeNames)
        if ~isempty(str)
            
            num = sscanf(str,'(%d, %d, %d)');
            dat.seedPosition(y-1,:) = num';
            
        end
    end
end


if exist('MPN','var')
save([MPN 'dat.mat'],'dat')
end














