% clear all buttom
    set(h.push_clear,'callback',{@Call_Back_clear,h});
    % menu popup
    set(h.popup_menu,'callback',{@Call_Back_menu,h.subframe});
    % mode
    set(h.radio_mode(1:end),'callback',{@Call_Back_mode,h});
    % check show3D buttom
    set(h.check_show3D,'callback',{@Call_Back_show3D,h});
%     % transparency slider
%     set(h.slider_show3D,'callback',{@Call_Back_slider_show3D,h});
    % grid number
    set(h.edit_domain.grid_num(1:3),'callback',{@Call_Back_xyz_grid,h});
    % lattice constant
    set(h.edit_domain.lattice_constant(1:6),'callback',{@Call_Back_lattice_constant,h});
    % shape parameters
    set(h.edit_domain.shape(1:2),'callback',{@Call_Back_shape_parameter,h});
    % shape popup
    set(h.popup_shape,'callback',{@Call_Back_shape,h});
    % permittivity and permeability
    set(h.edit_material(1),'callback',{@Call_Back_permittivity,h});
    set(h.edit_material(2),'callback',{@Call_Back_permeability,h});
    set(h.edit_material(3),'callback',{@Call_Back_reciprocity,h});
    set(h.edit_material(4),'callback',{@Call_Back_chirality,h});
    % path popups
    set(h.edit_lattice_show,'callback',{@Call_Back_lattice,h});
    % load, save, run, close
    set(h.push(1),'callback',{@Call_Back_run,h});
    set(h.push(2),'callback',{@Call_Back_load,h});
    set(h.push(3),'callback',{@Call_Back_close,h});
    set(h.push(4),'callback',{@Call_Back_save,h})
    % eigenmode
    set(h.hpush_eigenmode,'callback',{@Call_Back_eigenmode,h});