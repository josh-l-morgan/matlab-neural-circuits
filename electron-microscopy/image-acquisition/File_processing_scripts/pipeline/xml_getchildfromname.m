function c=xml_getchildfromname(data,name)
%returns the first immediate child from data with the given name, if found
%By Daniel Berger for MIT-BCS Seung, July 13 2009

c=[];
for child=1:1:size(data.children,2)
  if strcmp(data.children(child).name,name)
    c=data.children(child);
    return;
  end;
end;