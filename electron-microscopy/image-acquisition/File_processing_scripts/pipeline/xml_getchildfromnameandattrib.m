function c=xml_getchildfromnameandattrib(data,name,attribnames,attribvalues)

c=[];
for child=1:1:size(data.children,2)
  if strcmp(data.children(child).name,name)
    %This node has the right name; check all attribute values
    go=1;
    for attrib=1:1:size(attribnames,2)
      if (strcmp(attribnames{attrib},data.children(child).attributes(attrib).name)==0) ...
        || (strcmp(attribvalues{attrib},data.children(child).attributes(attrib).value)==0) 
        go=0;
      end;
    end;
    if go==1
      c=data.children(child);
      return;
    end;
  end;
end;