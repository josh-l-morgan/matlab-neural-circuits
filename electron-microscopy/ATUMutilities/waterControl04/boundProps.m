function[wP] = boundProps(handles)


global watProps

wP = watProps;

if wP.y1<1
    wP.y1 = 1;
end

if wP.x1<1
    wP.x1 = 1;
end

if wP.y2>wP.iHeight
    wP.y2 = wP.iHeight;
end

if wP.x2>wP.iWidth
    wP.x2 = wP.iWidth;
end

if wP.win2<= wP.win1;
    wP.win2 = watProps.win1+.001;
    set(handles.slider_win2,'Value',wP.win2)
end