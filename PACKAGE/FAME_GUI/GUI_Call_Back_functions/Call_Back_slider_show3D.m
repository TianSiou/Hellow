function Call_Back_slider_show3D(hobj,event,h)
    
    c_shape      = get( h.popup_shape,'value');
    transparency = get( h.slider_show3D, 'value' );
    lattice_constant = get( h.edit_domain.lattice_constant, 'string' );
    lattice_constant = str2num(lattice_constant);
    
    axes(h.ax_show3D)
    cla(h.ax_show3D)
    
    switch c_shape 
        case 1 % SC
            sphere_radius   = str2num(get(h.edit_domain.shape(1),'string'));
            cylinder_radius = str2num(get(h.edit_domain.shape(2),'string'));
            Plot_B_SC(lattice_constant,sphere_radius,cylinder_radius,transparency,h.ax_show3D);
        case 2 % SW
            width = str2num(get(h.edit_domain.shape(4),'string'));
            Plot_B_semiwoodpile(lattice_constant, width, transparency, h.ax_show3D);
        case 3 % FCC
            sphere_radius   = str2num(get(h.edit_domain.shape(6),'string'));
            cylinder_radius = str2num(get(h.edit_domain.shape(7),'string'));
            Plot_B_FCC(lattice_constant,sphere_radius,cylinder_radius,transparency,h.ax_show3D);
        case 4 % BCC
            iso_num1   = str2num(get(h.edit_domain.shape(8),'string'));
            iso_num2   = str2num(get(h.edit_domain.shape(9),'string'));
            Plot_B_BCC(lattice_constant,iso_num1,iso_num2,transparency,h.ax_show3D);    
    end
end