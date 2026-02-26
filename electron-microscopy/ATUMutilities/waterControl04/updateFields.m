function[] = updateFields(handles)

global watProps

set(handles.axes_Vid,'xtick',[],'ytick',[]);
set(handles.axes_Plot,'xtick',[],'ytick',[]);
set(handles.text_sourceDirectory,'string',watProps.sourceDir);
set(handles.edit_boxHeight,'string',watProps.boxHeight)
set(handles.edit_boxWidth,'string',watProps.boxWidth)
set(handles.edit_integrationTime,'string',watProps.intTime)
set(handles.slider_win1,'Value',watProps.win1)
set(handles.slider_win2,'Value',watProps.win2)
set(handles.slider_upDown,'Value',1-watProps.y1/watProps.iHeight)
set(handles.slider_leftRight,'Value',watProps.x1/watProps.iWidth)
set(handles.slider_thresh1,'Value',watProps.thresh1/256);
set(handles.edit_thresh1,'string',watProps.thresh1)
set(handles.edit_pumpInterval,'string',watProps.pumpInterval)
set(handles.edit_pumpDuration,'string',watProps.pumpDuration)
set(handles.togglebutton_threshPump,'Value',watProps.autoOn)
set(handles.togglebutton_sound,'Value',watProps.soundOn)
set(handles.listbox_colorType,'Value',watProps.colorMode)