function FAME_Tools_Display_Info( Pmesh, Plattice, Pmaterial, Par_recip_lattice_wave_vec_array )
    fprintf('----------------Computational Information----------------\n')
    fprintf('\tGrid numbers: [%d %d %d]\n',Pmesh.grid_num);
    fprintf(['\tLattice Type: ', Plattice.lattice_type, '\n']);
    fprintf('\tLattice Constants:\n ');
    fprintf('\t\t\ta     = %f\n',Plattice.lattice_constant_orig.a)
    fprintf('\t\t\tb     = %f\n',Plattice.lattice_constant_orig.b)
    fprintf('\t\t\tc     = %f\n',Plattice.lattice_constant_orig.c)
    fprintf('\t\t\talpha = %f\n',Plattice.lattice_constant_orig.alpha)
    fprintf('\t\t\tbeta  = %f\n',Plattice.lattice_constant_orig.beta)
    fprintf('\t\t\tgamma = %f\n',Plattice.lattice_constant_orig.gamma)
    fprintf(['\tMaterial Type: ', Pmaterial.material_type, '\n']);
    switch Pmaterial.material_type
        case 'isotropic'
            fprintf('\tPermittivity = [ ')
            for i = 1:length(Pmaterial.ele_permitt_in)
                fprintf('%.2f ',Pmaterial.ele_permitt_in(i))
                fprintf(']\n');
            end
            fprintf('\tPermeability = [ ')
            for i = 1:length(Pmaterial.mag_permeab_in)
                fprintf('%.2f ',Pmaterial.mag_permeab_in(i))
                fprintf(']\n');
            end
        case 'anisotropic'  
            for i = 1:length(Pmaterial.ele_permitt_in)
                fprintf('\tPermittivity(%d) = \n',i)
                display(Pmaterial.ele_permitt_in{i})
            end
            for i = 1:length(Pmaterial.mag_permeab_in)
                fprintf('\tPermeability(%d) = \n',i)
                display(Pmaterial.mag_permeab_in{i})
            end
        case 'biisotropic'
            fprintf('\tPermittivity = [ ')
            for i = 1:length(Pmaterial.ele_permitt_in)
                fprintf('%.2f ',Pmaterial.ele_permitt_in(i))
                fprintf(']\n');
            end
            fprintf('\tPermeability = [ ')
            for i = 1:length(Pmaterial.mag_permeab_in)
                fprintf('%.2f ',Pmaterial.mag_permeab_in(i))
                fprintf(']\n');
            end
            fprintf('\tReciprocity  = [ ')
            for i = 1:length(Pmaterial.mag_permeab_in)
                fprintf('%.2f ',Pmaterial.reciprocity_in(i))
                fprintf(']\n');
            end
            fprintf('\tChirality    = [ ')
            for i = 1:length(Pmaterial.mag_permeab_in)
                fprintf('%.2f ',Pmaterial.chirality_in(i))
                fprintf(']\n');
            end
    end
%     fprintf('\tBrillouin zone path: %s\n',path_string(1:end-1))
%     fprintf('\tPath partition number: %d\n',Plattice.part_num)
    fprintf('\tWave vector number: %d\n',size(Par_recip_lattice_wave_vec_array,2))
    fprintf('------------------------------------------------------\n')
end