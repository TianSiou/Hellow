function Call_Back_clear(hobj,event,h)
    axes(h.ax_display_shape_primitive)
    cla reset
    axes(h.ax_display_shape_unit)
    cla reset
    set(h.list_current_material_sphere, 'string', {});
    set(h.list_current_material_cylinder, 'string', {});
    set(h.list_parameter,   'string', {});
    set(h.element_parameter,'string',{});
    set(h.element_color,    'string',{});
    set(h.element_permitt,      'string',{});
    set(h.element_permeab,      'string',{});
    set(h.element_xi,      'string',{});
    set(h.element_zeta,      'string',{})
    set(h.element_edit{1},      'string',{});
    set(h.element_edit{2},      'string',{});
    set(h.element_edit{3},      'string',{});
    set(h.element_edit{4},      'string',{});
    set(h.element_edit{5},      'string',{});
    set(h.element_edit{6},      'string',{});
    set(h.list_parameter_material, 'string' ,{});
    set(h.object_tempinfo_sphere_center, 'string' ,{});
    set(h.object_tempinfo_sphere_material, 'string' ,{});
    set(h.object_tempinfo_sphere_radius, 'string' ,{});
    set(h.object_tempinfo_cylinder_material, 'string' ,{});
    set(h.object_tempinfo_cylinder_bot, 'string' ,{});
    set(h.object_tempinfo_cylinder_top, 'string' ,{});
    set(h.object_tempinfo_cylinder_radius, 'string' ,{});
    
end