function[colStruc] = col2struc(fileName);


%% read in col file
 fid = fopen(fileName);
  tline = fgetl(fid);
  y=1;
  while ischar(tline)
    if ((numel(tline)>0)&&(tline(1)~='%'))
      disp(tline);
      [a,count,errmsg,nextindex]=sscanf(tline, '%d   %d    %d %d %d %d   %d %d %d %d   %d %d %d   %d %d %d %d   %d   %d %d %d %d %d %d   ');
      data(y,:)=a';
      n=tline(nextindex:end); %Rest of line is name
      n=n(2:end-1); %Remove "" from name
      name{y}=n;
      y=y+1;
    end;
    tline = fgetl(fid);
  end;
  %mydata = textscan(fid, '%d   %d    %d %d %d %d   %d %d %d %d   %d %d %d   %s');
  fclose(fid);
  
  firstelement=data(1,17);
  data=data(2:end,:);
  name=name(2:end);

  %%
  colStruc.names = name;
  colStruc.data = data;
  colStruc.id = data(:,1);
  colStruc.seed = data(:,11:13);
  
  

