function Call_Back_shape(hobj,event,h)
% �� shape �U�Ԧ���檺 call back function
    bcolor             = [1 1 1 ]/1.06;
    hframe_shapeoption = h.sssframe;
    hframe_show3D      = h.frame_show3D;
    
    c_shape = get(h.popup_shape,'value'); % c_shape = 1,2,3 ���O�� SC, SW, FCC
    
%% �� Domain ��檺�ܰ�
    % structure parameter ���ܰ�
    for i = 1:4 %����������, ����A���}��ܤ� shpae �ҹ�����
        set( hframe_shapeoption(i),'visible','off');
    end
    set( hframe_shapeoption(c_shape),'visible','on'); % ���}��ܤ� shape �ҹ����� option frame
    
    hax = get(hframe_show3D,'child'); 
    cla(hax)
    if get(h.check_show3D, 'value')
        set(hax,'cameraviewangle',10);
    end
    
%     x_edge_len   = str2num(get(h.edit_domain(4),'string'));
    lattice_constant = get( h.edit_domain.lattice_constant, 'string' );
    lattice_constant = str2num(lattice_constant);
    
    transparency = get( h.slider_show3D, 'value' );
    % �� frame(1) �Pshow 3D ���ܰ�
    switch c_shape
        case 1
            ttext           = ' (Simple Cubic)';
            sphere_radius   = str2num(get(h.edit_domain.shape(1),'string'));
            cylinder_radius = str2num(get(h.edit_domain.shape(2),'string'));
            Plot_B_SC(lattice_constant, sphere_radius, cylinder_radius, transparency, hax);
            view(-39,31)
            
            set(h.edit_domain.grid_num(1:3),'string',num2str(16));
        case 2
            ttext = ' (Semi Woodpile)';
            width = str2num(get(h.edit_domain.shape(4),'string'));
            Plot_B_semiwoodpile(lattice_constant, width,transparency,hax);
            view(-45,25)
            
            set(h.edit_domain.grid_num(1:3),'string',num2str(16));
        case 3
            ttext           = ' (FCC)'; 
            sphere_radius   = str2num(get(h.edit_domain.shape(6),'string'));
            cylinder_radius = str2num(get(h.edit_domain.shape(7),'string'));
%             a_fcc           = sqrt(2)*x_edge_len;
            Plot_B_FCC(lattice_constant,sphere_radius,cylinder_radius,transparency,hax);
            view(-23,19)
            
            set(h.edit_domain.grid_num(1:3),'string',num2str(12));
        case 4
            ttext           = ' (BCC)'; 
            iso_num1   = str2num(get(h.edit_domain.shape(8),'string'));
            iso_num2   = str2num(get(h.edit_domain.shape(9),'string'));

            view(-23,19)
            Plot_B_BCC(lattice_constant,iso_num1,iso_num2,transparency,h.ax_show3D);    
            set(h.edit_domain.grid_num(1:3),'string',num2str(12));    
    end
    c_mode = cell2mat( get(h.radio_mode(1:end),'value') );
    c_mode = find(c_mode == 1);
    ttext = [ get(h.radio_mode(c_mode),'string'), ttext ];

    set(h.log_text(1),'string',ttext); % ��ܥثe�ҿ�ܪ� mode �P 3D structure
%% �� lattice ��檺�ܰ�
    
    switch c_shape
        case {1, 2}
            set( h.edit_lattice_show, 'string', 'GXMRG' )
            set(h.text_lattice_SC1(1:4),'visible','on');
            set(h.text_lattice_SC2(1:4),'visible','on');
            set(h.text_lattice_FCC1(1:6),'visible','off');
            set(h.text_lattice_FCC2(1:6),'visible','off');
            set(h.text_lattice_BCC1(1:7),'visible','off');
            set(h.text_lattice_BCC2(1:7),'visible','off');
            
        case 3
            set( h.edit_lattice_show, 'string', 'XULGWKX' )
            set(h.text_lattice_SC1(1:4),'visible','off');
            set(h.text_lattice_SC2(1:4),'visible','off');
            set(h.text_lattice_FCC1(1:6),'visible','on');
            set(h.text_lattice_FCC2(1:6),'visible','on');
            set(h.text_lattice_BCC1(1:7),'visible','off');
            set(h.text_lattice_BCC2(1:7),'visible','off');
            
        case 4
            set( h.edit_lattice_show, 'string', 'LHPNZQJL' )
            set(h.text_lattice_SC1(1:4),'visible','off');
            set(h.text_lattice_SC2(1:4),'visible','off');
            set(h.text_lattice_FCC1(1:6),'visible','off');
            set(h.text_lattice_FCC2(1:6),'visible','off');
            set(h.text_lattice_BCC1(1:7),'visible','on');
            set(h.text_lattice_BCC2(1:7),'visible','on');
            
    end
end