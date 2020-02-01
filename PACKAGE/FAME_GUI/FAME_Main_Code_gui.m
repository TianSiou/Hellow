function Pgraph = FAME_Main_Code_gui(par)

fprintf('--------------------------------Start FAME program----------------------------------\n')
fprintf('Loading package folder path......')
FAME_Folder_Manager( )
fprintf('Done\n')
%% Setting parameters
fprintf('Loading user options......')
[ Popt ] = FAME_User_Option_gui(par);
fprintf('Done\n')
%% Generate lattice vectors
grid_nums = [Popt.mesh.x_grid_num Popt.mesh.y_grid_num Popt.mesh.z_grid_num];
[ Popt.lattice.lattice_vec_a, Popt.mesh.edge_lens, Popt.mesh.mesh_lens, Omega, Popt.lattice.lattice_type, Popt.lattice.lattice_constant, Popt.material ] = ...
    FAME_Parameter_Generator( grid_nums, Popt.lattice.lattice_type, Popt.lattice.lattice_constant, Popt.material );
%% Partition the path of Brillouin zone 
fprintf('Generate Brillouin zone path point......')
[ wave_vec_array, wave_vec_num, path_string, path_string_new ] = ...
    FAME_Parameter_Brillouin_Zone_Path( Popt.lattice.part_num, Popt.lattice.lattice_type, Popt.lattice.lattice_constant, Popt.lattice.lattice_vec_a );
fprintf('Done\n')
%% (Main Step 1) - Construct the matrices B's
display_material = 'off';
fprintf('Construct material dependent matrix B......')
switch Popt.material.material_type
    case 'isotropic'
    [ B_eps, B_mu, invB_eps, invB_mu ] = ...
        FAME_Matrix_B_Isotropic( Popt.material.ele_permitt, Popt.material.mag_permeab, grid_nums, Popt.mesh.mesh_lens, Popt.lattice.lattice_vec_a, Popt.lattice.lattice_constant, Omega, Popt.material.material_handle, display_material, []);
    case 'anisotropic'
    [ B_eps, B_mu, invB_eps, invB_mu ] = ...
        FAME_Matrix_B_Anisotropic( Popt.material.ele_permitt, Popt.material.mag_permeab, grid_nums, Popt.mesh.mesh_lens, Popt.lattice.lattice_vec_a, Popt.lattice.lattice_constant, Omega, Popt.material.material_handle, display_material, [] );
    case 'biisotropic'
    [  B_eps,  B_mu,  B_xi,  B_zeta  ] = ...
        FAME_Matrix_B_Biisotropic( Popt.material.ele_permitt, Popt.material.mag_permeab, Popt.material.reciprocity, Popt.material.chirality, grid_nums, Popt.mesh.mesh_lens, Popt.lattice.lattice_vec_a, Popt.lattice.lattice_constant, Omega, Popt.material.material_handle, display_material, [] );    
end
fprintf('Done\n')
% Set waitbar data
hax_waitbar = par.hax_waitbar;
axes(hax_waitbar);
rotate3d off
cla reset;
xline = [0 1 1 0]; yline = [0 0 1 1];
l = line(xline,yline,'Color','k','EraseMode','none'); 
h = patch([0,1,1,0],[0,0,1,1],'r','EdgeColor','k','EraseMode','none');  
set(hax_waitbar,'XTick',[],'YTick',[]);
% Âk¹s
xpatch = [0 0 0 0];
set(h,'XData',xpatch);
title('Please wait... 0%');
drawnow

% Setting for Pgraph
node_num = length(path_string_new);
disconti_node_idx = [];
if node_num > 1
    for path_idx = 2 : node_num 
        node_idx = (path_idx-1)*(Popt.lattice.part_num-1) + 1;
        if length(path_string_new{path_idx}) == 3
            disconti_node_idx = [disconti_node_idx, node_idx];
        end
    end
end
xdata = 0;

% disconti_node_idx
% setdiff(1:wave_vec_num-1,disconti_node_idx)

for i = setdiff(1:wave_vec_num-1,disconti_node_idx)
    xdata(i+1) = norm(wave_vec_array(:,i+1) - wave_vec_array(:,i));
end
xdata = cumsum(xdata)
xtick = xdata(1:Popt.lattice.part_num-1:end)

xticklabal = path_string_new

% Display compuational information
FAME_Tools_Display_Info(grid_nums,Popt.lattice,path_string,Popt.material,wave_vec_num)
fprintf('===========================Start computing band structure===========================\n')
for i = 1:wave_vec_num
    % Choose wave vector in Brillouin zone path
    wave_vec = wave_vec_array(:,i);
    %% (Main Step 2) - Construct Lambdas corresponding to the choosen wave vector 
    fprintf('Construct Lambdas...')
    Lambdas = FAME_Matrix_Lambdas( wave_vec, grid_nums, Popt.mesh.mesh_lens, Popt.lattice.lattice_type, Popt.lattice.lattice_constant, Popt.lattice.lattice_vec_a );
    fprintf('Done\n')
    %% (Main Step 3) - Use fast algorithms to compute the smallest frequency for given wave vector, material type, and lattice type
    fprintf('Start FAME(%s, wave vector = [%.2f,%.2f,%.2f])...',Popt.material.material_type,wave_vec)
    switch Popt.material.material_type
        case 'isotropic'
        switch Popt.lattice.lattice_type
            case 'simple_cubic'
            [ Freq_array(:,i), Ele_field_mtx(:,:,i), cpu_time(i), LS_iter{i}, LS_cpu_time{i} ] = ...
                FAME_Fast_Algorithms_Simple_Isotropic(Popt.mesh.x_grid_num, Popt.mesh.y_grid_num, Popt.mesh.z_grid_num, wave_vec, B_eps, Lambdas, Popt.graph.plot_num);
            otherwise
            [ Freq_array(:,i), Ele_field_mtx(:,:,i), cpu_time(i), LS_iter{i}, LS_cpu_time{i}] = ...
                FAME_Fast_Algorithms_General_Isotropic(Popt.mesh.x_grid_num, Popt.mesh.y_grid_num, Popt.mesh.z_grid_num, wave_vec, B_eps, Lambdas, Popt.graph.plot_num);
        end
        [ C , Cs ] = FAME_Matrix_Curl( wave_vec, grid_nums, Popt.mesh.edge_lens, Popt.mesh.mesh_lens, {Popt.boundary.bd_cond_x,Popt.boundary.bd_cond_y,Popt.boundary.bd_cond_z}, Popt.lattice.lattice_type, Popt.lattice.lattice_vec_a, Popt.lattice.lattice_constant );
        Err(:,i) = FAME_Error_Check_Isotropic( Freq_array(:,i), Ele_field_mtx(:,:,i),C, Cs, B_eps );
        case 'anisotropic'
        switch Popt.lattice.lattice_type
            case 'simple_cubic'
            [ Freq_array(:,i), Ele_field_mtx(:,:,i), cpu_time(i), LS_iter{i}, LS_cpu_time{i} ] = ...
                FAME_Fast_Algorithms_Simple_Anisotropic(Popt.mesh.x_grid_num, Popt.mesh.y_grid_num, Popt.mesh.z_grid_num, wave_vec, invB_eps, Lambdas, Popt.graph.plot_num);
            otherwise
            [ Freq_array(:,i), Ele_field_mtx(:,:,i), cpu_time(i), LS_iter{i}, LS_cpu_time{i} ] = ...
                FAME_Fast_Algorithms_General_Anisotropic(Popt.mesh.x_grid_num, Popt.mesh.y_grid_num, Popt.mesh.z_grid_num, wave_vec, invB_eps, Lambdas, Popt.graph.plot_num);
        end
        [ C , Cs ] = FAME_Matrix_Curl( wave_vec, grid_nums, Popt.mesh.edge_lens, Popt.mesh.mesh_lens, {Popt.boundary.bd_cond_x,Popt.boundary.bd_cond_y,Popt.boundary.bd_cond_z}, Popt.lattice.lattice_type, Popt.lattice.lattice_vec_a, Popt.lattice.lattice_constant );
        Err(:,i) = FAME_Error_Check_Anisotropic( Freq_array(:,i), Ele_field_mtx(:,:,i),C, Cs, B_eps );
        case 'biisotropic'
        switch Popt.lattice.lattice_type
            case 'simple_cubic'
            [ Freq_array(:,i), EleMag_field_mtx(:,:,i), cpu_time(i), LS_iter{i}, LS_cpu_time{i} ] = ...
                FAME_Fast_Algorithms_Simple_Biisotropic(Popt.mesh.x_grid_num, Popt.mesh.y_grid_num, Popt.mesh.z_grid_num, wave_vec, B_eps, B_mu, B_xi, B_zeta, Lambdas, Popt.graph.plot_num   );
            otherwise
            [ Freq_array(:,i), EleMag_field_mtx(:,:,i), cpu_time(i), LS_iter{i}, LS_cpu_time{i} ] = ...
                FAME_Fast_Algorithms_General_Biisotropic(Popt.mesh.x_grid_num, Popt.mesh.y_grid_num, Popt.mesh.z_grid_num, wave_vec, B_eps, B_mu, B_xi, B_zeta, Lambdas, Popt.graph.plot_num   );
        end
        [ C , Cs ] = FAME_Matrix_Curl( wave_vec, grid_nums, Popt.mesh.edge_lens, Popt.mesh.mesh_lens, {Popt.boundary.bd_cond_x,Popt.boundary.bd_cond_y,Popt.boundary.bd_cond_z}, Popt.lattice.lattice_type, Popt.lattice.lattice_vec_a, Popt.lattice.lattice_constant );
        Err(:,i) = FAME_Error_Check_Biisotropic( Freq_array(:,i), EleMag_field_mtx(:,:,i),C, Cs, B_eps, B_mu, B_xi, B_zeta );
    end
    fprintf('Done! Use time %.2f(sec.)\n',cpu_time(i))
    % Set waitbar data
    value = i/wave_vec_num;
    xpatch = [0 value value 0];
    set(h,'XData',xpatch);
    title(['Please wait... ',num2str(floor(100*value)),'%']);
    drawnow
    
    clear('Lambdas','wave_vec','C','Cs');
    fprintf('================================== %3.1f%% complete =====================================\n',100*(i/wave_vec_num))
end

title('Complete!')
axes(par.hax_band)
%% Plot band structure
% FAME_Plot_Band_Structure( path_string_new, Popt.lattice.part_num, Freq_array )
%% Collecting the data for plotting the band structure and eigen modes
% The order of nodes on wave vector lattice
%node_array = (0 : Plattice.node_num - 1) * (Popt.lattice.part_num - 1) + 1;
% Collecting the data for plotting the band structure and eigen modes
Pgraph.eigenvec_array = Ele_field_mtx;
Pgraph.eigvalue_array = Freq_array;
%% Collecting the data for plotting the band structure and eigen modes
% The order of nodes on wave vector lattice
% node_array = (0 : Plattice.node_num - 1) * (Popt.lattice.part_num - 1) + 1;
% Collecting the data for plotting the band structure and eigen modes
% Pgraph.eigenvec_array = eigenvec_array( : , 1 : Popt.graph.plot_num , : );
% Pgraph.eigvalue_array = eigvalue_array(     1 : Popt.graph.plot_num , : );
FAME_GUI_Graphic_Plotter_Band( Freq_array , path_string_new, Popt.lattice.part_num, xdata,xtick,xticklabal)
Pgraph.xdata      = xdata;
Pgraph.xtick      = xtick;
Pgraph.xticklabal = xticklabal;

Pgraph.wave_vec_array = wave_vec_array;
Pgraph.sn             = [Popt.mesh.x_grid_num Popt.mesh.y_grid_num Popt.mesh.z_grid_num 3];
Pgraph.par            = par;
Pgraph.path_string    = path_string_new;
Pgraph.part_num       = Popt.lattice.part_num;
% save savedata Popt Pgraph par
