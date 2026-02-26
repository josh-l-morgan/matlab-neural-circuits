function ResGui
%Display and edit results

%% Make Figure
RG = figure('Visible','off','Position',[10,10,1100,520]);
Display = axes('Units','pixels','Position',[10,10,500,500]);

%% Folder Managment
GetFolder = uicontrol('Style','pushbutton',...
    'String','Get Folder',...
    'Callback',{@GetFolder_Callback},...
    'Position',[520,490,80,20]);
ShowFile= uicontrol('Style','text','Position',[600,490,220,20],...
    'String','No File Selected','FontSize',12);
Save= uicontrol('Style','pushbutton','Position',[820,490,70,20],...
    'String','Save','Callback',{@Save_Callback});



%% Manage Cell Data
CellFrame = uicontrol('Style','frame','Position',[520,430,370,50]);
CellFrameStr = uicontrol('Style','text','Position',[530,445,60,30],...
    'String','Cell Info','FontSize',10);
AgeStr = uicontrol('Style','text','Position',[650,440,50,15],...
    'String','Cell Age');
Age= uicontrol('Style','edit','Position',[705,440,40,15],...
    'String','?');
TypeStr = uicontrol('Style','text','Position',[650,460,50,15],...
    'String','Cell Type');
Type= uicontrol('Style','edit','Position',[705,460,40,15],...
    'String','?');
ON = uicontrol('Style','pushbutton','Position',[750,460,30,15],...
    'String','ON','Callback',{@ON_Callback});
OFF = uicontrol('Style','pushbutton','Position',[785,460,30,15],...
    'String','OFF','Callback',{@OFF_Callback});
BI = uicontrol('Style','pushbutton','Position',[820,460,30,15],...
    'String','BI','Callback',{@BI_Callback});
Other = uicontrol('Style','pushbutton','Position',[855,460,30,15],...
    'String','Other','Callback',{@Other_Callback});

ClassStr = uicontrol('Style','text','Position',[765,440,50,15],...
    'String','Cell Class');
Class= uicontrol('Style','edit','Position',[815,440,60,15],...
    'String','?');

ReadyStr = uicontrol('Style','text','Position',[530,440,50,15],...
    'String','Ready?');
Ready= uicontrol('Style','edit','Position',[580,440,15,15],...
    'String','?');
ReadyYes = uicontrol('Style','pushbutton','Position',[600,440,20,15],...
    'String','Y','Callback',{@ReadyYes_Callback});
ReadyNo = uicontrol('Style','pushbutton','Position',[620,440,20,15],...
    'String','N','Callback',{@ReadyNo_Callback});

Notes= uicontrol('Style','edit','Position',[520,10,380,50]);
NotesStr= uicontrol('Style','text','Position',[530,40,30,10],...
    'String','Notes:');


%% Display Results

Rx=520;Ry=320; %% set lower corner for group
StatsFrame = uicontrol('Style','frame','Position',[Rx,Ry,240,100]);
StatsHead = uicontrol('Style','text','Position',[Rx+60,Ry+75,160,20],...
    'String','  Length       Dots       Dots/Dend');
StatsCell = uicontrol('style','text','Position',[Rx+5,Ry+55,50,20],...
    'String','Cell');
StatsA1 = uicontrol('style','text','Position',[Rx+5,Ry+30,50,20],...
    'String','Arbor1');
StatsA2 = uicontrol('style','text','Position',[Rx+5,Ry+5,50,20],...
    'String','Arbor2');


%%Run analysis
CBx=910; CBy=470;
FindThreshCB = uicontrol('Style','checkbox','Position',[CBx,CBy,80,15],...
    'String','FindThresh');
SkeletonizeCB = uicontrol('Style','checkbox','Position',[CBx,CBy-15,80,15],...
    'String','Skeletonize');
DotFindCB = uicontrol('Style','checkbox','Position',[CBx,CBy-30,80,15],...
    'String','DotFind');
RatioCB = uicontrol('Style','checkbox','Position',[CBx,CBy-45,80,15],...
    'String','Ratio');
MaskCB = uicontrol('Style','checkbox','Position',[CBx,CBy-60,80,15],...
    'String','Mask');
DrawCB = uicontrol('Style','checkbox','Position',[CBx,CBy-75,80,15],...
    'String','Draw');
DotDendCB = uicontrol('Style','checkbox','Position',[CBx,CBy-90,80,15],...
    'String','DotDend');
ShollCB = uicontrol('Style','checkbox','Position',[CBx,CBy-105,80,15],...
    'String','Sholl');
FixDataCB = uicontrol('Style','checkbox','Position',[CBx,CBy-120,80,15],...
    'String','FixData');
Run = uicontrol('Style','pushbutton','Position',[CBx,CBy-150,80,25],...
    'String','Run','FontSize',12,...
    'Callback',{@Run_Callback});

%%Map




%% Initialize the GUI.
% Change units to normalized so componets resize automatcally.
set([GetFolder,ShowFile,Display,ON,OFF,BI,Other,Notes,...
    Save,AgeStr,TypeStr,Age,Type],'Units','normalized');


set(RG,'Name','Result Viewer')
movegui(RG,'center')
set(RG,'Visible','on') %make Gui visible
if ~exist('.\temp'),mkdir('.\temp'),end

%%Map Raid
Raid='\\128.208.64.36\wonglab\Josh\Analyzed\Output';
RaidDir=dir(Raid);
RD=struct2cell(RaidDir);
RD=RD(1,:); RD=RD(3:size(RD,2));
SelectCellRaid = uicontrol('Style','text','Position',[1010,455,80,30],...
    'String','Select Cell From Raid');
SelectCell = uicontrol('Style','listbox','Position',[1010,10,80,450],...
    'String', RD,...
    'Callback',{@SelectCell_Callback});  
  


%% Initialize Variables (check previous)
global Status Results Cell TPN DPN

if exist('.\temp\RG_LastFile.mat')
    load(['.\temp\RG_LastFile.mat'])
    Back=find(SPN=='\');
    DFN=SPN(Back(size(Back,2))+1:size(SPN,2))
    TPN=[SPN '\'];
    DPN=FindImage(TPN);
    [Cell,Results,Status,Combo]=LoadData(SPN);
    Cell,Results,Status,
    set(RG,'CurrentAxes',Display)
    image(Combo)
    ShowResults
end


%% Callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Folder management
    function  GetFolder_Callback(source,eventdata)
        Save_Callback
        clear Cell Results Status
        axes(Display),image(0);
        Clear
        %Get Folder (default to previous)
        if exist('SPN') , if ischar(SPN),    SPN=uigetdir(SPN,'*.*');
        else,    SPN=uigetdir('*.*'); end
        else,     SPN=uigetdir('*.*');   
        end
        save('.\temp\RG_LastFile.mat','SPN')
        Back=find(SPN=='\');
        DFN=SPN(Back(size(Back,2))+1:size(SPN,2))
        TPN=[SPN '\'];
        DPN=FindImage(TPN);        
        [Cell,Results,Status,Combo]=LoadData(SPN);
        Cell,Results,Status
        set(RG,'CurrentAxes',Display)
        image(Combo)
        ShowResults
    end

    function Save_Callback(source,eventdata)
        if exist('DFN')
            Cell.Name=DFN;
            Cell.Notes=get(Notes,'String');
            Cell.Age=get(Age,'String');
            Cell.Type=get(Type,'String');
            Cell.Ready=get(Ready,'String');
            Cell.Class=get(Class,'String');
            save([SPN '\Cell.mat'],'Cell');
            Cell
        end
    end

%% Cell data Management
    function ON_Callback(source,eventdata)
        set(Type,'String','ON');end
    function OFF_Callback(source,eventdata)
        set(Type,'String','OFF');end
    function BI_Callback(source,eventdata)
        set(Type,'String','BI');end
    function Other_Callback(source,eventdata)
        set(Type,'String','Other');end
    %%Set Ready State
    function ReadyYes_Callback(source,eventdata)
        set(Ready,'String','Y');end
    function ReadyNo_Callback(source,eventdata)
        set(Ready,'String','N');end

 
%% Display Results
    function ShowResults(source,eventdata)
        DPN
         ShowFile= uicontrol('Style','text','Position',[600,490,200,20],...
            'String',DFN,'FontSize',12);
        
        %%Set Available Variables
        if isstruct(Cell)
            set(Type,'String',Cell.Type);
            set(Age,'String',Cell.Age);
            set(Ready,'String',Cell.Ready);
            set(Notes,'String',Cell.Notes);
            set(Class,'String',Cell.Class);
        else
            set(Type,'String','?');
            set(Age,'String','?');
            set(Ready,'String','?');  
            set(Notes,'String','?');
            set(Class,'String','?');
        end
        
        %%Check Results
        if isstruct(Results)
            if isfield(Results.CellStats,'TotalLength')
                CLength=Results.CellStats.TotalLength;
                CDots=Results.CellStats.TotalDots;
                CDD=round(1000*Results.CellStats.AverageDensity)/1000;
                A1Length=Results.Arbor(1).Length;
                A1Dots=Results.Arbor(1).Dots;
                A1DD=round(1000*Results.Arbor(1).DotDend)/1000;
            
            
               if size(Results.Arbor,2)>1 
                   A2Length=Results.Arbor(2).Length;
                   A2Dots=Results.Arbor(2).Dots;
                   A2DD=round(1000*Results.Arbor(2).DotDend)/1000;
               else
                    A2Length=0;
                   A2Dots=0;
                   A2DD=0;
               end
               
            else
                CLength=0;  CDots=0;  CDD=0;   A1Length=0;  A1Dots=0;
                A1DD=0;   A2Length=0;   A2Dots=0;  A2DD=0;    NN=0;
                L=0;   DDxy=0;   DDdepth=0;
             end % if cell stats total length
            
           % if neares Neighbor has been run
           if isfield(Results.CellStats,'NN')
               NN=Results.CellStats.NN;
           else,    NN=0;      end
           
           % if Depth has been run
           if isfield(Results,'Depth')
               L=Results.Depth.DendBinDepth; 
               DDdepth=Results.Depth.DotPerDendDepth; 
           else
               L=0;
               DDdepth=0;
           end
           
           %If XY was run
           if isfield(Results,'XY')
               DDxy=Results.XY.DotDend; 
           else
               DDxy=0;
           end

        else  %If no results exist
            CLength=0;  CDots=0;  CDD=0;   A1Length=0;  A1Dots=0;
            A1DD=0;   A2Length=0;   A2Dots=0;  A2DD=0;    NN=0;
            L=0;   DDxy=0;   DDdepth=0;
        end  %Finished checking results
        
           LengthCell = uicontrol('style','text','Position',[Rx+70,Ry+55,40,20],...
                 'String',round(CLength));
           DotCell = uicontrol('style','text','Position',[Rx+120,Ry+55,40,20],...
                 'String',CDots);
           DotDendCell = uicontrol('style','text','Position',[Rx+170,Ry+55,40,20],...
                'String',CDD);
            
           LengthA1 = uicontrol('style','text','Position',[Rx+70,Ry+30,40,20],...
                'String',round(A1Length));
           DotA1 = uicontrol('style','text','Position',[Rx+120,Ry+30,40,20],...
                'String',A1Dots);
           DotDendA1 = uicontrol('style','text','Position',[Rx+170,Ry+30,40,20],...
                'String',A1DD);
           
          LengthA2 = uicontrol('style','text','Position',[Rx+70,Ry+5,40,20],...
                  'String',round(A2Length));
          DotA2 = uicontrol('style','text','Position',[Rx+120,Ry+5,40,20],...
                  'String',A2Dots);
          DendA2 = uicontrol('style','text','Position',[Rx+170,Ry+5,40,20],...
                  'String',A2DD);        
          
           
           %%Plot Nearest Neighbor

               NNh=hist(NN);
               NNplot = axes('Unit','pixels','Position',[Rx+265,Ry+20,100,70]);
               NNplotStr = uicontrol('Style','text','Position',[Rx+270,Ry+90,90,15],...
                   'String','Nearest Neighbor');
               bar(NNplot,NNh);
               LessThen2=round((sum(NN<=2)/size(NN,2))*100);
               NNLT2 = uicontrol('Style','text','Position',[Rx+310,Ry+50,25,15],...
                   'String',LessThen2);
           
           
           %%Plot Dend Depth

               Lplot = axes('Unit','pixels','Position',[690,255,200,60]);
               plot(Lplot,L)
               LplotStr = uicontrol('Style','text','Position',[540,275,100,35],...
                   'FontSize',10,'String','Dendritic Length Across Depth');

               DDdepthPlot = axes('Unit','pixels','Position',[690,170,200,60]);
               plot(DDdepthPlot,DDdepth)          
               DDdepthPlotStr = uicontrol('Style','text','Position',[540,190,100,35],...
                   'FontSize',10,'String','Dot Density Across Depth');  
           
           

               DDxyPlot = axes('Unit','pixels','Position',[690,85,200,60]);
               plot(DDxyPlot,DDxy)        
               DDxyPlotStr = uicontrol('Style','text','Position',[540,95,110,50],...
                'FontSize',10,'String','Dot Density with distance from Cell Body');       

        %% Display status
        StaFrame = uicontrol('Style','frame','Position',[910,10,80,220]);
        FieldNames=fieldnames(Status);
        FNsize=size(FieldNames,1);
        for i=1:FNsize;
            y=40+(180/FNsize)*(FNsize-i);
            StatName = uicontrol('Style','text','Position',[915,y,40,15],...
                'String',FieldNames(i));
            ShowStat = uicontrol('Style','text','Position',[965,y,15,15],...
                'String',getfield(Status,char(FieldNames(i))));
            
        end
        
        
    end %% End Show Results 




    %% Check Box Callbacks
    %{
    function FindThresh_Callback(hObject, eventdata, handles)
    function Skeletonize_Callback(hObject, eventdata, handles)
    function Ratio_Callback(hObject, eventdata, handles)
    function Mask_Callback(hObject, eventdata, handles)
    function Draw_Callback(hObject, eventdata, handles)
    function DotDend_Callback(hObject, eventdata, handles)
    function Sholl_Callback(hObject, eventdata, handles)
    function FixData_Callback(hObject, eventdata, handles)
    function DotFind_Callback(hObject, eventdata, handles)
    %}

    function Run_Callback(source,eventdata)
        Save_Callback      
        hrtext = uicontrol('Style','text','String','Running','Position',[910,300,80,20]);
        pause(.01)
         Status.Running=1; ShowResults; pause(.01)
         plot(Display,0)
        if get(FindThreshCB,'Value'); anaFT; end
        if get(SkeletonizeCB,'Value'); anaSk; end
        if get(DotFindCB,'Value'); anaDF; end
        if get(FixDataCB,'Value'); anaFix;end
        if get(RatioCB,'Value'); anaRa; end
        if get(MaskCB,'Value'); anaMa; end
        if get(DrawCB,'Value'); anaDr; end
        if get(DotDendCB,'Value'); anaDD(TPN); end
        if get(ShollCB,'Value'); anaSh(TPN); end
        if exist([TPN 'data\DotStats.mat']), anaNN; end
        Status.Running=0;
       
        [Cell,Results,Status,Combo]=LoadData(SPN);
        Save_Callback;
        ShowResults
        Display = axes('Units','pixels','Position',[10,10,500,500]);
        image(Combo)
        
        hrtext = uicontrol('Style','text','String','Done','Position',[910,300,80,20]);
        pause(.01)
       
    end

    function SelectCell_Callback(source,eventdata)
        Save_Callback
        clear Cell Results Status
        axes(Display),image(0);
        TargCell=RD(get(SelectCell,'Value'))
        SPN=[Raid '\' char(TargCell)]
        save('.\temp\RG_LastFile.mat','SPN')
        Back=find(SPN=='\');
        DFN=SPN(Back(size(Back,2))+1:size(SPN,2))
        [Cell,Results,Status,Combo]=LoadData(SPN);
        Cell,Results,Status
        set(RG,'CurrentAxes',Display)
        image(Combo)
        ShowResults
        TPN=[SPN '\'];
        DPN=FindImage(TPN);
    
    end
   

    function Clear(source,eventdata)
        'Clear'
    end
    

end