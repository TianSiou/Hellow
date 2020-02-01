function Call_Back_menu(hobj,~,hframe) 
%給下拉式選單(1)所使用的 call back function
    for i = 1:length(hframe)
        set( hframe(i),'visible','off','position',[.01 .01 .01 .01])
    end
    i = get(hobj,'value');
    set( hframe(i),'visible','on','position',[.05 .15 .9 .83])
end