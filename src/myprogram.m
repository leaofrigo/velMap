function myprogram

    f = figure('WindowStyle','normal');
    ax = axes;
    x = 0:100;
    y = x.^2;

    plotline = plot(x,y);
    c = uicontextmenu;

    % Set c to be the plot line's UIContextMenu
    plotline.UIContextMenu = c;

    % Create menu items for the uicontextmenu
    m1 = uimenu(c,'Label','dashed','Callback',@setlinestyle);
    m2 = uimenu(c,'Label','dotted','Callback',@setlinestyle);
    m3 = uimenu(c,'Label','solid','Callback',@setlinestyle);

        function setlinestyle(source,callbackdata)
            switch source.Label
                case 'dashed'
                    plotline.LineStyle = '--';
                case 'dotted'
                    plotline.LineStyle = ':';
                case 'solid'
                    plotline.LineStyle = '-';
            end
        end
end
