function lattice_constant = FAME_Parameter_Lattice_Constants_Format( varargin )
    lattice_type = varargin{1};
    if isstruct(varargin{2}) == 1
        lattice_constant = varargin{2};
        lattice_constant = check_lattice_constant( lattice_type, lattice_constant );
    else
        lattice_vec_a = varargin{2};
        lattice_constant = check_lattice_vec( lattice_type, lattice_vec_a );
    end

end

function lattice_constant = check_lattice_constant( lattice_type, lattice_constant )

% This routine reture pretreatment lattice constants
    lattice_constant.old = lattice_constant;
    switch lattice_type  
        case {'simple_cubic','face_centered_cubic','body_centered_cubic'}
              %% Cubic system
        % In cubic system, lattice constant must satisfying
        %               a = b = c, £\=£]=£^=£k/2
            if isfield(lattice_constant,{'a','b','c'}) ~= 1
                error('There has only one lattice constant ''a'' should be specified. Notice that, in cubic system, lattice constant must satisfying a = b = c, £\=£]=£^=£k/2');
            end

            lattice_constant.b     = lattice_constant.a;
            lattice_constant.c     = lattice_constant.a;
            lattice_constant.alpha = pi/2;
            lattice_constant.beta  = pi/2;
            lattice_constant.gamma = pi/2;
        case {'hexagonal' ,'rhombohedral'}
              %% Hexagonal and Rhombohedral systems
        % In hexagonal and rhombohedral systems, lattice constant must satisfying
        %               a = b, £\=£]=£k/2, £^=2£k/3
            if sum(isfield(lattice_constant,{'a','c'})) ~= 2
                error('The lattice constant ''a'' or ''c'' were not specified. Notice that, in hexagonal and rhombohedral systems, lattice constant must satisfying a = b, £\=£]=£k/2, £^=2£k/3');
            end
            lattice_constant.b     = lattice_constant.a;
            lattice_constant.alpha = pi/2;
            lattice_constant.beta  = pi/2;
            lattice_constant.gamma = 2*pi/3;
        case {'primitive_tetragonal','body_centered_tetragonal'}
              %% Tetragonal system
        % In tetragonal system, lattice constant must satisfying
        %               a = b ¡Ú c, £\=£]=£^=£k/2
            if sum(isfield(lattice_constant,{'a','c'})) ~= 2
                error('The lattice constant ''a'' or ''c'' were not specified. Notice that, in tetragonal system, lattice constant must satisfying a = b ¡Ú c, £\=£]=£^=£k/2');
            elseif lattice_constant.a == lattice_constant.c
                error('The lattice constants ''a'' and ''c'' must be distinct. Notice that, in tetragonal system, lattice constant must satisfying a = b ¡Ú c, £\=£]=£^=£k/2');
            end
            lattice_constant.b     = lattice_constant.a;
            lattice_constant.alpha = pi/2;
            lattice_constant.beta  = pi/2;
            lattice_constant.gamma = pi/2;
        case {'primitive_orthorhombic','a_base_centered_orthorhombic','c_base_centered_orthorhombic','face_centered_orthorhombic', 'body_centered_orthorhombic'}
              %% Orthorhombic system
        % In orthorhombic system, lattice constant must satisfying
        %               a ¡Ú b ¡Ú c, £\=£]=£^=£k/2    
            if sum(isfield(lattice_constant,{'a','b','c'})) ~= 3
                error('The lattice constant ''a'', ''b'' or ''c'' were not specified. Notice that, in orthorhombic system, lattice constant must satisfying a ¡Ú b ¡Ú c, £\=£]=£^=£k/2');
            elseif (lattice_constant.a == lattice_constant.b) || (lattice_constant.a == lattice_constant.c) || (lattice_constant.c == lattice_constant.b)
                error('The lattice constants ''a'', ''b'' and ''c'' must be distinct. Notice that, in orthorhombic system, lattice constant must satisfying a ¡Ú b ¡Ú c, £\=£]=£^=£k/2');
            end
            lattice_constant.alpha = pi/2;
            lattice_constant.beta  = pi/2;
            lattice_constant.gamma = pi/2;
        case {'primitive_monoclinic','base_centered_monoclinic'}
              %% Monoclinic system
        % In monoclinic system, lattice constant must satisfying
        %               a, b <= c, £\<£k/2, £]=£^=£k/2        
            if isfield(lattice_constant,{'alpha'}) ~= 1
                error('There are only £\ were be specified. Notice that, in monoclinic system, lattice constant must satisfying a, b <= c, £\<£k/2, £]=£^=£k/2         ')
            elseif lattice_constant.alpha >= pi/2
                error('The lattice constant £\ were be smaller than £k/2. Notice that, in monoclinic system, lattice constant must satisfying a, b <= c, £\<£k/2, £]=£^=£k/2         ')
            elseif (lattice_constant.a > lattice_constant.c) || (lattice_constant.b > lattice_constant.c)
                error('The lattice constant a,b,c must satisfied a, b ? c. Notice that, in monoclinic system, lattice constant must satisfying a, b <= c, £\<£k/2, £]=£^=£k/2         ')
            end

            lattice_constant.beta  = pi/2;
            lattice_constant.gamma = pi/2;
        case 'triclinic'    
             %% Triclinic system
       % In triclinic system, lattice constant must contains all 6
       % constants
            if sum(isfield(lattice_constant,{'a','b','c','alpha','beta','gamma'})) ~= 6
                error('The lattice constant ''a'', ''b'', ''c'', ''£\'',''£]'' or ''£]'' were not specified.');
            end
    end
end

function lattice_constant = check_lattice_vec( lattice_type, lattice_vec_a )
    
    a1 = lattice_vec_a(:,1); a2 = lattice_vec_a(:,2); a3 = lattice_vec_a(:,3);

    switch lattice_type  
        case {'triclinic'}
            lattice_constant.a = norm(a1);
            lattice_constant.b = norm(a2);
            lattice_constant.c = norm(a3);
            lattice_constant.alpha = dot(a2,a3)/(lattice_constant.b*lattice_constant.c);
            lattice_constant.beta  = dot(a1,a3)/(lattice_constant.a*lattice_constant.c);
            lattice_constant.gamma = dot(a1,a2)/(lattice_constant.a*lattice_constant.b);
        otherwise
            error('This lattice type is not valid in current version!');
    end
end