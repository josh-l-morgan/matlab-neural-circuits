function param=photostagemap_ovregionabsolutealign(param)
%This function computes an absolute coordinates in a global coordinate
%space for all overview slice images
%By Daniel Berger for MIT-BCS Seung / Harvard Lichtman, March 2010

lparam=param.ovregionabsolutealign;

for w=param.processwafers
  %Compute which slices have been updated manually on this wafer by comparing original and partially corrected stagemaps
  originalstagemap=[param.stagemapsdir sprintf(lparam.originalstagemap,w)];
  modifiedstagemap=[param.stagemapsdir sprintf(lparam.partiallycorrectedstagemap,w)];
  differentslices=find(comparestagemaps(originalstagemap,modifiedstagemap));
  
  nrofslices=size(param.alignoverviewimages.stsec{w},1); %nrofslices=size(param.generatestagemap.linearslicepos{w},1);
  gavail=param.ovregionalign.gavail{w};
  gtransx=param.ovregionalign.gtransx{w};
  gtransy=param.ovregionalign.gtransy{w};
  grot=param.ovregionalign.grot{w};
  gcorr=param.ovregionalign.gcorr{w};
  gmatrix=param.ovregionalign.gmatrix{w};
  
  gabsmatrix=cell(nrofslices,lparam.nrofpredictions);
  gabsmatrixcount=zeros(nrofslices,1);
  %gabsmatrix{1}=[[1 0 0]; [0 1 0]; [0 0 1]]; gabsmatrixcount(1)=1;
  for i=1:1:size(differentslices,1)
    gabsmatrix{differentslices(i)}=[[1 0 0]; [0 1 0]; [0 0 1]]; gabsmatrixcount(differentslices(i))=1;
  end;
  
  %Do recursive alignment to slice 1
  switch lparam.alignmode
    case 0  % lparam.npredictions predictions for each slice, average all equally
      %while sum(gabsmatrixcount==0)>0
      while sum(gabsmatrixcount==lparam.nrofpredictions)<nrofslices
        %for rec=1:1:nrofrecursions
        for i=1:1:nrofslices
          if gabsmatrixcount(i)>0 %source exists
            for j=1:1:nrofslices
              if gavail(i,j) %&&(gcorr(i,j)>lparam.mincorr)) %transformation source-target exists and correlation is high enough
                if gabsmatrixcount(j)<lparam.nrofpredictions
                  gabsmatrixcount(j)=gabsmatrixcount(j)+1;
                  %gabsmatrix{j,gabsmatrixcount(j)}=gabsmatrix{i,1}*squeeze(gmatrix(i,j,:,:));
                  gabsmatrix{j,gabsmatrixcount(j)}=squeeze(gmatrix(i,j,:,:))*gabsmatrix{i,1};
                end;
              end;
            end;
          end;
        end;
      end;
      %average predictions
      for i=1:1:nrofslices
        ggabsmatrix{i}=gabsmatrix{i,1};
        if lparam.nrofpredictions>1
          %compute arithmetic mean
          for j=2:1:lparam.nrofpredictions
            ggabsmatrix{i}=ggabsmatrix{i}+gabsmatrix{i,j};
          end;
          ggabsmatrix{i}=ggabsmatrix{i}/lparam.nrofpredictions;
          ggabsmatrix{i}=makerigidmatrix(ggabsmatrix{i}); %normalize to make rigid transformation
        end;
      end;
  end;
  param.ovregionabsolutealign.gabsmatrixcount{w}=gabsmatrixcount;
  param.ovregionabsolutealign.gabsmatrix{w}=ggabsmatrix;
end;