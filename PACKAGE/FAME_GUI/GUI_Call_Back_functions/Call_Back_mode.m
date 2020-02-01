function Call_Back_mode(hobj,~,h)
    
    c_mode         = cell2mat( get(h.radio_mode(1:end),'value'));
    c_mode         = find(c_mode == 1);
    material_type  = get(h.radio_mode(c_mode),'string');

    FAME_GUI_Info_Set('material_type', material_type, h)
    [Popt] = FAME_GUI_Info_Get(h);
    
    
    i = get(h.popup_shape,'value');
    %% 對狀態顯示的變動
    set(h.log_text(1),'string',[Popt.material.material_type,' (',Popt.lattice.lattice_type,')']); % 顯示目前所選擇的 mode 與 3D structure  
    
    %% 對 material 選單的變動
    switch c_mode 
        case 1
            set(h.text_material(3:4),'visible','off')
            set(h.edit_material(3:4),'visible','off')
            set(h.edit_material(1),'string','13');
            set(h.edit_material(2),'string','1');
        case 2
            set(h.text_material(3:4),'visible','on')
            set(h.edit_material(3:4),'visible','on')
            set(h.edit_material(1),'string','13');
            set(h.edit_material(2),'string','1');
        case 3
            set(h.text_material(3:4),'visible','off')
            set(h.edit_material(3:4),'visible','off')
            set(h.edit_material(1),'string','[21.26 0 -14i; 0 16 0; 14i 0 21.26]');
            set(h.edit_material(2),'string','[1 0 0; 0 1 0; 0 0 1]');
    end
end