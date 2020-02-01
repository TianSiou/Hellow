function Call_Back_show3D(hobj,event,h)
% µ¹ show 3D structure «öÁäªº call back function
    i = get(hobj,'value');
    h.isShow3D = i;
    display_unit = 0;
    display_prim = 0;
    if i == 1
        if display_unit == 1
            set( h.frame_show3D_unit,'visible','on');
            set( h.frame_show3D_prim,'visible','off');
        else
            set( h.frame_show3D_unit,'visible','off');
            set( h.frame_show3D_prim,'visible','on');
        end
    else
        set( h.frame_show3D_unit,'visible','off');
        set( h.frame_show3D_prim,'visible','off');
    end
end