function tentative2
sld = uicontrol('Style', 'slider',...
        'Min',1,'Max',50,'Value',41,...
        'Position', [400 20 120 20],...
        'Callback', @surfzlim); 
    function surfzlim(source,callbackdata)
%         val = 51 - source.Value;
        % For R2014a and earlier:
        val = get(source,'Value')
        get(source)


    end
end