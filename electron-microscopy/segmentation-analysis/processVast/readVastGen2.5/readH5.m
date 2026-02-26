


SPN = 'D:\LGNs1\Segmentation\rhoanna\run3combo\'
MPN = [SPN(1:end-1) '_mat']
if ~exist(MPN), mkdir(MPN),end

%"/label_index/123"

%% explore h5
h5name = dir([SPN '*.h5']);
h5file = [SPN h5name.name];
%h5disp(h5File)
h5map = h5info(h5file)

%%

numOb = length(h5info.Groups.Datasets);

testSet = h5readatt(h5file,['/Datasets/' h5info.Groups.Datasets(1).Name]);




for i = 1:length(useIds)
    o = useIds(i);
    vastSubs{o} = zeros(idSize(i),3,'uint16');
end





%save([TPN 'vastSubs.mat'],'vastSubs','-v7.3')

