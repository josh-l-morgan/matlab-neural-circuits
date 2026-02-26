function[vastSubs] = fixVastSubs(vastSubs)

global glob


%%Embarassing hak
if strcmp(glob.NA.export.exportName,'Josh_B') & ...
        strcmp(glob.vol.activeName,'v1')
    hackCellNavGlob
end


if isfield(glob.NA.export,'fix')

    %%Shift z position of all subs
    if isfield(glob.NA.export.fix,'changeZ')
        disp('Shifting z for vastSubs')
        for i = 1:length(vastSubs)
            if ~isempty(vastSubs{i})
                vastSubs{i}(:,3) = vastSubs{i}(:,3) + glob.NA.export.fix.changeZ;
            end
        end
    end

end

