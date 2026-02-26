function [Files,labels]=GetSortedImageList(directory)

filestruct = dir([directory filesep '*.tif']);
%labels=zeros(1,length(filestruct));
%Files=cell(1,length(filestruct));
%MatFiles=cell(1,length(filestruct));
  soNum = 0;
for i = 1:length(filestruct)
            %Extract Label
            nam  = filestruct(i).name;
            listNum = [];
            for s = 1:length(nam)
                foundNum = str2num(nam(s));
                if foundNum > .5
                    listNum  = [listNum  nam(s)];
                end
            end
                
            Label = listNum;
          
            if  str2num(Label)>0
                soNum = soNum+1;
                Files{soNum}=[directory filesep filestruct(i).name];
                %MatFiles{soNum}=[directory filesep filestruct(i).name(1:end-3) 'mat'];
                labels(soNum) = str2num(Label);
            end
end
[labels,indices]=sort(labels);
Files=Files(indices);
%MatFiles=MatFiles(indices);