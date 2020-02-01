function Call_Back_run(hobj,event,h)
% 抓取參數並執行FAME main code
    global Pgraph

    gpar.hfig             = h.fig;
    gpar.hax_waitbar = h.ax_waitbar;
    gpar.hax_band    = h.ax_band;
%% Multi-core 設定
%     gpar.isMultiCore = get( h.check_muticore,'value');

%% mode 設定 
    c_mode = cell2mat( get(h.radio_mode(1:end),'value'));
    c_mode = find(c_mode == 1);
    switch c_mode
        case 1 
            gpar.mode = 'isotropic';
        case 2 
            gpar.mode = 'biisotropic';
        case 3 
            gpar.mode = 'anisotropic';
    end
%% Domain 參數設定
    %Domain 設定
    gpar.x_grid_num = str2num(get(h.edit_domain.grid_num(1),'string'));
    gpar.y_grid_num = str2num(get(h.edit_domain.grid_num(2),'string'));
    gpar.z_grid_num = str2num(get(h.edit_domain.grid_num(3),'string'));
    gpar.lattice_constant.a     = str2num(get(h.edit_domain.lattice_constant(1),'string'));
    gpar.lattice_constant.b     = str2num(get(h.edit_domain.lattice_constant(2),'string'));
    gpar.lattice_constant.c     = str2num(get(h.edit_domain.lattice_constant(3),'string'));
    gpar.lattice_constant.alpha = str2num(get(h.edit_domain.lattice_constant(4),'string'));
    gpar.lattice_constant.beta  = str2num(get(h.edit_domain.lattice_constant(5),'string'));
    gpar.lattice_constant.gamma = str2num(get(h.edit_domain.lattice_constant(6),'string'));
    %3D structure parameter 設定
    gpar.sphere_radius    = str2num(get(h.edit_domain.shape(1),'string'));
    gpar.cylinder_radius  = str2num(get(h.edit_domain.shape(2),'string'));

    lattice_type = get(h.log_text(1),'string');
    temp_idx1 = find(lattice_type == '(');
    temp_idx2 = find(lattice_type == ')');
    gpar.lattice_type = lattice_type(temp_idx1+1:temp_idx2-1);
%% material 參數設定
    popup_num                     = get(h.popup_shape,'value');
    file_name_list                = get(h.popup_shape,'string');
    gpar.file_name                = file_name_list(popup_num,:);
    gpar.file_name(isspace(gpar.file_name)) = [];
    gpar.ele_permitt = str2num( get(h.edit_material(1),'string') );
    gpar.mag_permeab = str2num( get(h.edit_material(2),'string') );
    gpar.reciprocity = str2num( get(h.edit_material(3),'string') );
    gpar.chirality   = str2num( get(h.edit_material(4),'string') );
%% lattice 路徑
%     cc = get(h.edit_lattice_show,'string');
%     gpar.path=lower(cc);
%     gpar.lattice_part_num = str2num( get(h.edit_lattice_part_num,'string') );
%% eigensolver 設定
    gpar.lattice_part_num = str2num( get(h.edit_lattice_part_num,'string') );
    c_eigensolver =  get(h.radio_eigensolver(1),'value');
    if c_eigensolver == 1
        gpar.eig_solver = 'jdsiraSEP';
    else
        gpar.eig_solver = 'eigs';
    end
%% run FAME main code
    %防止 band structure 顯示在 show 3D 上
    set(h.frame_show3D_unit,'visible','off');%關閉 show 3D structure 的frame
    set(h.frame_show3D_prim,'visible','off');
    set(h.check_show3D,'value',0);
    Call_Back_show3D( h.check_show3D,[],h );
%     set(h.ax_band,'visible','on');%打開顯示band structure 的frame
    cla(h.ax_band);axes(h.ax_band);

    Pgraph = FAME_Main_Code_gui(gpar);

%% 標出 band gap      
    minx = min(Pgraph.eigvalue_array');
    maxx = max(Pgraph.eigvalue_array');
    nn   = max(Pgraph.eigvalue_array(:));
    xx   = [0:nn/1000:nn];
    for i = 1:size(Pgraph.eigvalue_array,1)
        xx(find(xx>=minx(i) & xx<= maxx(i)))=0;
    end
    im  = ones(length(xx),20);
    idx = find(xx~=0);

    if ~isempty(idx)
        im(idx,:) = 0;
        figure(1)
        xL = xlim;
        im(:,:,2) = 1;
        im(:,:,3) = im(:,:,1);
        im(:,:,1) = 1;
        hold on
        xxx = [0:nn/1000:nn];
        cla;
        image([xL(1): (xL(2)-xL(1))/19:xL(2)],xxx,im);

        bgap = xx(idx(end))-xx(idx(1));
        text(0.5,(xx(idx(1))+xx(idx(end)))/2,[ 'Band Gap=' num2str(bgap,'%6.3g') ],'fontsize',15,'color','r','FontWeight','bold');
        wv = Pgraph.xdata;%[[0:12] 12+[1:4]*(3*.5^2)^.5/.5];
        plot(wv,Pgraph.eigvalue_array','-bo','Linewidth',2)
        node_num = length(Pgraph.path_string);
        if node_num > 1
            for path_idx = 2 : node_num 
                node_idx = (path_idx-1)*(Pgraph.part_num-1) + 1;
                if length(Pgraph.path_string{path_idx}) == 3
                    plot([wv(node_idx),wv(node_idx)],[0,max(max(Pgraph.eigvalue_array))],'k--','linewidth', 2)
                end
            end
        end
        alpha(0.5);
    end

%% 打開 eigenmode 的按鍵

    set(h.hpush_eigenmode,'visible','on');
end