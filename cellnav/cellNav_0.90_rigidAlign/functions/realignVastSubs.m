function[vastSubs] = realignVastSubs(vastSubs)

global glob


load("MPN.mat");
if exist([MPN 'shiftZ.mat'],'file')
    disp('shifting planes')
    load([MPN 'shiftZ.mat']);
        for v = 1:length(vastSubs)
            disp(sprintf('fixing sub %d of %d',v,length(vastSubs)))
            if ~isempty(vastSubs{v})
                vastSubs{v} = rigidPtsStack(vastSubs{v},shiftZ.As);
            end
        end
end




