function[dat] = parseGoogleDat()


load('MPN.mat')
load([MPN 'dat.mat'])

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


for x = 1:size(rawDat,2);
   
        str = rawDat{1,x};
        if isempty(str)
            str = sprintf('col %d',x);
        end
        dat.columnNames{x} = str;
end

dat.rawDat = rawDat;

%% Parse IDs
for y = 2:size(rawDat,1)
   
    str = rawDat{y,1};
    nums = sscanf(str,'%d');
    if ~isempty(nums)
        dat.cid(y-1) = nums(1);
    end

end


%% parse cell type
dat.typeNames = {'rgc' 'tcr' 'lin' 'unk','phot','hrz','bpc','amc','mul'};
for y = 2:size(rawDat,1)
    
   dat.type(y-1) = 0; 
   str = lower(rawDat{y,2});
   for t = 1:length(dat.typeNames)
       if sum(regexp(str,dat.typeNames{t}))
         dat.type(y-1) = t;
         break
       end
   end
           
end


%% parse on/off
dat.onOffNames = {'on' 'off'};
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
dat.subTypeNames{8} = {'vgc'};
dat.subTypeNames{7} = {'bc1' 'bc2' 'bc3a' 'bc3b' 'bc4' 'bc5i' 'bc5o' 'bc5t' 'bc6' 'bc7' 'bc8' 'bc11' 'unk'};

for y = 2:size(rawDat,1)
    
    str = lower(rawDat{y,4});
    if isempty(str)
        dat.subType(y-1) = 0;
    else
        checkNames = dat.subTypeNames{dat.type(y-1)};
        for t = 1:length(checkNames)
            if sum(strcmp(str,checkNames{t}))
                dat.subType(y-1) = t;
                break
            end
        end
    end
end




%% parse tracing
dat.tracingNames = {'none' 'fragment' 'partial' 'full' 'error'};
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



save([MPN 'dat.mat'],'dat')














