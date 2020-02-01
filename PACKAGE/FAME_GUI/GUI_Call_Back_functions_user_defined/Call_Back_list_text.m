function Call_Back_list_text(hobj,event,h)
    answer = inputdlg({'Parameter name:','Value:'},'Add parameter');
    if isempty(answer) ~= 1
        par_list = get(h.list_parameter,'string');
        par_list{end+1} = [answer{1} ,'=',answer{2},';'];
        set(h.list_parameter,'string',par_list);
        set(h.list_parameter,'value',length(par_list));
        par_list_material = get(h.list_parameter_material,'string');
        par_list_material{end+1} = ['material.scalar.',answer{1},'=',answer{2},';'];
        set(h.list_parameter_material,'string',par_list_material);
       
    end
end