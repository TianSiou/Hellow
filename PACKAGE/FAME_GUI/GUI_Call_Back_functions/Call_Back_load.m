function Call_Back_load(hobmatj,event,h)
    global Pgraph

    [filename,pathname] = uigetfile('*.mat');
    if filename ~= 0
        try
            load([pathname filename]);
            try
                fprintf(['start load file ', filename,'\n'])
                
                set(h.ax_band,'visible','on');%ゴ秨陪ボband structure frame
                set(h.hpush_eigenmode,'visible','on');

                % 砞﹚ mode
                if strcmp(par.mode,'isotropic') == 1
                    set(h.radio_mode(1),'value',1);
                else
                    set(h.radio_mode(2),'value',1);
                end
                % 砞﹚ domain 把计
                set( h.edit_domain.grid_num(1), 'string', num2str(par.x_grid_num) );
                set( h.edit_domain.grid_num(2), 'string', num2str(par.y_grid_num) );
                set( h.edit_domain.grid_num(3), 'string', num2str(par.z_grid_num) );
                
                set( h.edit_domain.shape(1), 'string', num2str(par.sc_sphere_radius) );
                set( h.edit_domain.shape(2), 'string', num2str(par.sc_cylinder_radius) );
                set( h.edit_domain.shape(3), 'string', num2str(par.sw_wood_length) );
                set( h.edit_domain.shape(4), 'string', num2str(par.sw_wood_width) );
                set( h.edit_domain.shape(5), 'string', num2str(par.sw_wood_height) );
                set( h.edit_domain.shape(6), 'string', num2str(par.fcc_sphere_radius) );
                set( h.edit_domain.shape(7), 'string', num2str(par.fcc_cylinder_radius) );
                if strcmp( par.shape, 'simple_cubic' ) == 1
                    set( h.popup_shape, 'value', 1 );
                    c_shape = 1;
                elseif strcmp( par.shape, 'semi_woodpile' ) == 1
                    set( h.popup_shape, 'value', 2 );
                    c_shape = 2;
                elseif strcmp( par.shape, 'face_centered_cubic' ) == 1
                    set( h.popup_shape, 'value', 3 );
                    c_shape = 3;
                end
                % 砞﹚ material 把计
                set( h.edit_material(1), 'string', num2str(par.ele_permitt) );
                set( h.edit_material(2), 'string', num2str(par.mag_permeab) );
                set( h.edit_material(3), 'string', num2str(par.reciprocity) );
                set( h.edit_material(4), 'string', num2str(par.chirality) );
                % 砞﹚ lattice path
                switch c_shape
                    case {1,2}
                        set( h.text_lattice_SC(1:4),'visible','on' )
                        set( h.text_lattice_FCC(1:6),'visible','off' )
                    case 3
                        set( h.text_lattice_SC(1:4),'visible','off' )
                        set( h.text_lattice_FCC(1:6),'visible','on' )
                end
                path = upper( par.path );
                text_path = [];
                for i = 1:length(path)-1
                    text_path = [text_path, path(i),'->'];
                end
                text_path = ['Path: ', text_path, path(end)];
                set( h.text_lattice_show, 'string', text_path ); 
                % 砞﹚ eigensolver
                if strcmp( par.eig_solver, 'jdsiraSEP' ) == 1
                    set( h.radio_eigensolver(1), 'value', 1);
                elseif strcmp( par.eig_solver, 'eigs' ) == 1
                    set( h.radio_eigensolver(2), 'value', 1);
                end
            end
            % 礶 band structure
            cla(h.ax_band);axes(h.ax_band);
            FAME_Graphic_Plotter2( Popt.domain, Popt.material, Popt.lattice, Popt.graph, Pgrid, Plattice.wave_vec_num, Plattice.node_num, Pmaterial, Pgraph, Pgraph.xdata,Pgraph.xtick,Pgraph.xticklabal);
        end
    end
end