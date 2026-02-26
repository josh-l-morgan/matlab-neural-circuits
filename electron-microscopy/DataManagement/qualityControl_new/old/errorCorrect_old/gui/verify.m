function[proceed] = verify()
%Display and edit results

%% Make Figure
RG = figure('Visible','off','Position',[100,100,300,100]);

%% Buttons

ReadyStr = uicontrol('Style','text','Position',[10,50,280,50],...
    'String','Proced to Retake?','FontSize',20);
ReadyYes = uicontrol('Style','pushbutton','Position',[10,10,130,30],...
    'String','Yes','Callback',{@ReadyYes_Callback},'FontSize',20);
ReadyNo = uicontrol('Style','pushbutton','Position',[160,10,130,30],...
    'String','No','Callback',{@ReadyNo_Callback},'FontSize',20);


%% Initialize the GUI.

set(RG,'Name','Result Viewer')
movegui(RG,'center')
set(RG,'Menubar','none')
set(RG,'Visible','on') %make Gui visible

waitfor(RG, 'Visible','off') 
proceed = get(gcf, 'userdata');
close(RG)
%% Cell data Management

    %%Set Ready State
    function[answer] = ReadyYes_Callback(source,eventdata)
        set(gcf, 'userdata', 1); 
        set(RG, 'Visible','off') 
        
    end
    function[answer] =  ReadyNo_Callback(source,eventdata)
        set(gcf, 'userdata',0); 
        set(RG, 'Visible','off') 
    end

end