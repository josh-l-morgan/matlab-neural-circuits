function[qual] = q2qual(TPN);

if ~exist('TPN','var')
    TPN = GetMyDir;
end

if exist([TPN 'q.mat'])
    load([TPN 'q.mat'])


if exist([TPN 'mif.mat']);
    load([TPN 'mif.mat']);
else
    getMif(TPN);
end



checklist = find(ones(length(q.checkedNams),1));
for w = 1:length(mif.w)
    for s = 1:length(mif.w(w).sec)
        for t = 1:length(mif.w(w).sec(s).tileNams)
            foldnam = mif.w(w).sec(s).tileFolders(t) ;
            filenam =  mif.w(w).sec(s).tileNams(t);
            for o = 1:length(checklist)
                chk = checklist(o);
                if strcmp(q.checkedFolders(chk),foldnam);
                    if strcmp(q.checkedNams(chk),filenam);
                        r = mif.w(w).sec(s).rowIds(t);
                        c = mif.w(w).sec(s).colIds(t);
                        qual.w(w).sec(s).tile(r,c) = q.tile(chk).quality;
                        checklist = checklist(checklist~=chk);
                        break
                    end
                end
            end
        end
    end
end

save([TPN 'qual.mat'],'qual')
else
    'No quality file found'
end

