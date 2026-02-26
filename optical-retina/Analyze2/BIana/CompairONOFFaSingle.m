%% Compare ON vs OFF
%%different from previous ON OFF in that the only cut off is the dividing
%%line between arbors. 



    clear ONOFF
    TPN=GetMyDir
    load([TPN 'Use.mat'])
    
    %Get Rotation if it exists
    if exist([TPN 'Rot.mat']), 
        eRot=1; load([TPN 'Rot.mat']),else
        eRot=0; Rot=[1 0 0; 0 1 0; 0 0 1];
    end;
    DPos=CoRot(Use.DPos,Rot);
    Mids=CoRot(Use.Mids,Rot);
    
    
    load([TPN 'data\ResultsComp.mat'])
    
    
    %%Find cell position
      
    load([TPN 'data\ResultsComp.mat'])
    Top=[]; Bottom=[];
    for a = 1:size(Results.Arbor,2)
        Top(a)=Results.Arbor(a).Top;
        Bottom(a)=Results.Arbor(a).Bottom;
    end
    
       
    Cent=CoRot(Use.Cent,Rot);
    CellZ=Cent(3);
    if abs(CellZ-mean(Top))>abs(CellZ-mean(Bottom));
        Flip=1; else Flip=0; 
    end
    Fr(k)=Flip;
    
    Middle=(Bottom(2)+Top(1))/2;
    
    
    Layer=zeros(size(Mids,1),1);
    Layer((Mids(:,3)<Middle(1)))=1;
    if size(Top,2)>1
        Layer((Mids(:,3)>Middle) )=2;
    end

    DotLayer=zeros(size(DPos,1),1);
    DotLayer((DPos(:,3)>=Top(1)) & (DPos(:,3)<Middle(1)))=1;
    DotLayer((DPos(:,3)>Middle) )=2;



    ONOFF.DotLayer=DotLayer;
    ONOFF.MidLayer=Layer;
    
    save([TPN 'ONOFFa.mat'],'ONOFF')
    












