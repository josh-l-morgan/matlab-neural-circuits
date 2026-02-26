%% Turn mosaics from waffer into viewable (centers and downsample) mosaics

%% Get Waffer folder information
wif = GetMyWafer;

%% read log file
lnam = wif.log{end};

F=fopen(lnam,'r')
scanLog = textscan(F, '%s','Delimiter',';')
logF = scanLog{1};
fclose(F);

%% search log file

row = regexp(logF,' r=');
col = regexp(logF,' c=');
mos = regexp(logF,'Starting mosaic of site ');  %"w01_st01_sec03" at 6:49:06 PM
sid=zeros(length(logF),3);
r = zeros(length(logF),1);
c = zeros(length(logF),1);
rcInd = {};
for i = 1:length(logF)
    s = logF{i};

     if ~isempty(mos{i})
        q = strfind(s,'"');
        und = strfind(s,'_s');
        sec(1) = str2num(s(q(1)+3:und(1)-1));
        sec(2) = str2num(s(und(1)+3:und(2)-1));
        sec(3) = str2num(s(und(2)+4:q(2)-1));
        sid(i,:) = sec;
        
     elseif ~isempty(row{i})
        colend = regexp(s(col{i}:end),':');
        r(i) = str2num(s(row{i}+3:col{i}-2));
        c(i) = str2num(s(col{i}+3:col{i}+colend(1)-2)); 
        sid(i,:) = sec;
        
     elseif exist('sec','var')
         sid(i,:) = sec;
     end
end
%% Record indexes in cell as {waffer, strip, section, row, column}
Inds = cell(max(sid(:,1)),max(sid(:,2)),max(sid(:,3)),max(r),max(c));
for fw = 1:max(sec(:,1))
    for ft = 1:max(sec(:,2));
        for fs = 1:max(sec(:,3));
            for fr = 0:max(r)
                if ~fr
                    MosInd{fw,ft,fs} = find((sid(:,1) == fw)...
                        & (sid(:,2) == ft) & (sid(:,3) == fs)...
                        & (r==0));
                else
                    for fc = 1:fc
                        Inds{fw,ft,fs,fr,fc} = find((sid(:,1) == fw)...
                            & (sid(:,2) == ft) & (sid(:,3) == fs)...
                            & (r==fr)&(c==fc));
                    end
                end
            end
        end
    end
end

logSub = logF(Inds{1,1,1,1,1})
matLog.logF=logF;
matLog.mosInd = MosInd;
matLog.tileInd = Inds;
save([wif.dir 'matLog.mat'],'matLog')
