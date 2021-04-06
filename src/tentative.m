function tentative()

h=figure('Deletefcn',@func)
c=uicontrol(h,'Style','pushbutton','String','Pause','Callback',@pausefunc);
cpos = get(c,'Position');
c1pos = cpos;
c1pos(1) = c1pos(1)+c1pos(3);
uicontrol(h,'Style','pushbutton','String','Stop','Position',c1pos,'Callback',@stopfunc);
c2pos= c1pos;
c2pos(1) = c2pos(1)+c2pos(3);
uicontrol(h,'Style','pushbutton','String','run','Position',c2pos,'Callback',@runfunc);

    function func(source,event)
        play=false
    end
    function pausefunc(source,event)
        play=false
        pause(1);    
    end
    function stopfunc(source,event)
        play=false
        i=1
        pause(1);    
    end
    function runfunc(source,event)
        play=true
        pause(1);    
    end
end