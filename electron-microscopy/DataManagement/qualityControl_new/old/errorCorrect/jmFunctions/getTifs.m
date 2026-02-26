function[inams] = getTifs(TPN);

dTPN = dir(TPN); dTPN = dTPN(3:end);
inams = {};
%% get tifs
for i = 1:length(dTPN);
    nam = dTPN(i).name;
    if length(nam)>3
        if strcmp(nam(end-3:end),'.tif')
            inams = [inams nam];
        end
    end
end

%% sort tifs
% if ~isempty(inams)
% for i = 1:length(inams)
%     nam = inams{i};
%     nnams(i) = length(inams);  % set to last by default
%     for p = 1:length(nam)-4;
%         snam = nam(end - 3 - p:end -4);
%         nnam = str2num(snam);
%         if isempty('nnam')
%             break
%         end
%         nnams(i) = nnam;
%     end
% end
nnams = inams;
[sorted idx] = sort(nnams);
inams = inams(idx);
end
