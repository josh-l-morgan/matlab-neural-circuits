function[] = showStatus()

global watProps

if watProps.status == 2
    set(handles.text_Status,'backgroundcolor',[1 0 0 ]);
    set(handles.text_Status,'string','Watching, but NOT pumping')
elseif watProps.status == 3
    set(handles.text_Status,'backgroundcolor',[0 1 0]);
    set(handles.text_Status,'string','Auto pump is active')
else
    set(handles.text_Status,'backgroundcolor',[.7 .7 .4]);
    set(handles.text_Status,'string','confused')
end