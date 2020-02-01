%% varify Lattice type
switch get(h.pop_type,'value')
            case 1
                lattice_type     = 'simple_cubic';
            case 2
                lattice_type     = 'face_centered_cubic';
            case 3
                lattice_type     = 'body_centered_cubic';
            case 4
                lattice_type     = 'hexagonal';
            case 5
                lattice_type     = 'rhombohedral';
            case 6
                lattice_type     = 'primitive_tetragonal';
            case 7
                lattice_type     = 'body_centered_tetragonal';
            case 8
                lattice_type     = 'primitive_orthorhombic';
            case 9
                lattice_type     = 'face_centered_orthorhombic';
            case 10
                lattice_type     = 'body_centered_orthorhombic';
            case 11
                lattice_type     = 'a_base_centered_orthorhombic';
            case 12
                lattice_type     = 'c_base_centered_orthorhombic';
            case 13
                lattice_type     = 'primitive_monoclinic';
            case 14
                lattice_type     = 'a_base_centered_monoclinic';
            case 15
                lattice_type     = 'triclinic';
        end
