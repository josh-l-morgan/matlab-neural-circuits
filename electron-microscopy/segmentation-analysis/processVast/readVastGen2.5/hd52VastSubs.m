function[] = hd52VastSubs(SPN);


%SPN = 'D:\LGNs1\Segmentation\rhoanna\run3combo\'
MPN = [SPN(1:end-1) '_mat\']
if ~exist(MPN), mkdir(MPN),end

%"/label_index/123"

%% explore h5
if exist([MPN 'h5map.mat'],'file')
    load([MPN 'h5map.mat'])
    
else
h5name = dir([SPN '*.h5']);
h5file = [SPN h5name.name];
%h5disp(h5File)
h5map = h5info(h5file)
save([MPN 'h5map.mat'],'h5map')
end

%%

numOb = length(h5map.Groups.Datasets);
for i = 1:numOb
    disp(sprintf('reading %d of %d',i,numOb))
    datName = ['/label_index/' h5map.Groups.Datasets(i).Name];
    vastSubs{i} = h5read(h5file,datName)';
end

save([MPN 'vastSubs.mat'],'vastSubs','-v7.3')




