function Call_Back_run(hobj,event,h)
% 抓取參數並執行FAME main code
    gpar.hfig        = h.fig;
    gpar.hax_waitbar = h.ax_waitbar;
    gpar.hax_band    = h.ax_band;
%% Load user options from exist file "FAME_User_Option.m"
    [ Popt ] = FAME_GUI_Info_Get(h);
    Popt.lattice.lattice_type = 'user_defined';
%% Generate modified lattice vectors and lattice constants for computing
    [ Par.mesh, Par.lattice, Par.recip_lattice, Par.material, Par.eig ] = ...
        FAME_Parameter_Generator( Popt.mesh.grid_num, Popt.lattice.lattice_type, Popt.lattice.lattice_constant, Popt.material);
%% Generate wave vector array by partition the path of Brillouin zone 
    [ Par.recip_lattice ] = FAME_Parameter_Brillouin_Zone_Path( Popt.recip_lattice.part_num, Par.lattice, Par.recip_lattice );
%% Locating indices for the material inside
    [ Par.material.B.ele_x_idx, Par.material.B.ele_y_idx, Par.material.B.ele_z_idx, Par.material.B.mag_x_idx, Par.material.B.mag_y_idx, Par.material.B.mag_z_idx, Par.material.B.org_idx ] = ...
        FAME_Material_Locate_Index( Par.mesh, Par.lattice, Par.material);
%% Start FAME
    %防止 band structure 顯示在 show 3D 上
    set(h.frame_show3D_unit,'visible','off');%關閉 show 3D structure 的frame
    set(h.frame_show3D_prim,'visible','off');
    set(h.check_show3D,'value',0);
    Call_Back_show3D( h.check_show3D,[],h );
%     set(h.ax_band,'visible','on');%打開顯示band structure 的frame
    cla(h.ax_band);axes(h.ax_band);
    
    dir_name = ['result_',date]; % folder name
    temp     = dir;
    i = 1;
    while exist(dir_name) ~= 0
        dir_name = ['result_',date,'(',num2str(i),')']; % add series number to folder name if it has been used
        i = i + 1;
    end
    mkdir(dir_name); % make folder
    mfile_name                  = [dir_name,'\','PC_para'];            % set m-file name for 'Par'
    txtfile_name_para           = [dir_name,'\','PC_para.txt'];        % set text-file name for 'PC_para.txt'
    txtfile_name_B              = [dir_name,'\','invB.txt'];           % set text-file name for 'invB.txt'
    txtfile_name_wave_vec_array = [dir_name,'\','wave_vec_array.txt']; % set text-file name for 'wave_vec_array.txt'
    save(mfile_name, 'Par'); % save m-file
    FAME_mfile2txt_future_version( mfile_name, txtfile_name_para, txtfile_name_B, txtfile_name_wave_vec_array); % save text-files
    
    FAME_option.discrete_method = 'Yee_scheme';
    [ omega_array, ele_field, mag_field, comput_info ] = FAME_Main_Code( Par.mesh, Par.lattice, Par.material, Par.eig, Par.recip_lattice.wave_vec_array, FAME_option );

%% Plot band structure    
    omega_array = omega_array/(2*pi);
    FAME_Plot_Band_Structure( Par.recip_lattice.path_string, Par.recip_lattice.part_num, omega_array )

%% 標出 band gap      
    node_num = length(Par.recip_lattice.path_string_new);
    disconti_node_idx = [];
    if node_num > 1
        for path_idx = 2 : node_num 
            node_idx = (path_idx-1)*(Par.recip_lattice.part_num-1) + 1;
            if length(Par.recip_lattice.path_string_new{path_idx}) == 3
                disconti_node_idx = [disconti_node_idx, node_idx];
            end
        end
    end

%     xdata = 0;
%     for i = setdiff(1:Par.recip_lattice.wave_vec_num-1,disconti_node_idx)
%         xdata(i+1) = norm(Par.recip_lattice.wave_vec_array(:,i+1) - Par.recip_lattice.wave_vec_array(:,i));
%     end
%     xdata = cumsum(xdata);
%     xtick = xdata(1:Par.recip_lattice.part_num-1:end);

    minx = min(omega_array');
    maxx = max(omega_array');
    nn   = max(omega_array(:));
    xx   = [0:nn/1000:nn];
    for i = 1:size(omega_array,1)
        xx(find(xx>=minx(i) & xx<= maxx(i)))=0;
    end
    im  = ones(length(xx),20);
    idx = find(xx~=0);

    if ~isempty(idx)
        im(idx,:) = 0;
%         figure(1)
        xL = xlim;
        im(:,:,2) = 1;
        im(:,:,3) = im(:,:,1);
        im(:,:,1) = 1;
        hold on
        xxx = [0:nn/1000:nn];
%         cla;
        image([xL(1): (xL(2)-xL(1))/19:xL(2)],xxx,im);

        bgap = xx(idx(end))-xx(idx(1));
        text(0.5,(xx(idx(1))+xx(idx(end)))/2,[ 'Band Gap=' num2str(bgap,'%6.3g') ],'fontsize',15,'color','r','FontWeight','bold');
%         wv = xdata;%[[0:12] 12+[1:4]*(3*.5^2)^.5/.5];
% size(wv)
% size(omega_array')
%         plot(wv,omega_array','-bo','Linewidth',2)
%         node_num = length(Par.recip_lattice.path_string_new);
%         if node_num > 1
%             for path_idx = 2 : node_num 
%                 node_idx = (path_idx-1)*(Par.recip_lattice.part_num-1) + 1;
%                 if length(Par.recip_lattice.path_string_new{path_idx}) == 3
%                     plot([wv(node_idx),wv(node_idx)],[0,max(max(omega_array))],'k--','linewidth', 2)
%                 end
%             end
%         end
        alpha(0.5);
    end    
    
%% 打開 eigenmode 的按鍵
%     set(h.hpush_eigenmode,'visible','on');
end