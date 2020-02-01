function [lattice_vec_a,lattice_vec_a_orig,length_a1,length_a2,length_a3,theta_1,theta_2,theta_3,Permutation] = FAME_Parameter_Lattice_Vector(lattice_type,lattice_constant,lattice_vec_a)
% This routine return the lattice vectors and constants corresponding to
% given lattice type and lattice constants

%theta_1 is the angle between a1 and a2
%theta_2 is the angle between a1 and a3
%theta_3 is the angle between a2 and a3
    if exist('lattice_vec_a') ~= 1
        [a1_orig,a2_orig,a3_orig] = Lattice_vector_generate(lattice_type, lattice_constant);
    else
        a1_orig = lattice_vec_a(:,1);
        a2_orig = lattice_vec_a(:,2);
        a3_orig = lattice_vec_a(:,3);
    end
    
    
    lattice_vec_a_orig = [a1_orig,a2_orig,a3_orig];
    lattice_vec_a      = lattice_vec_a_orig;
    Permutation        = [1,2,3];
    
    length_a1_orig = norm(a1_orig);
    length_a2_orig = norm(a2_orig);
    length_a3_orig = norm(a3_orig);

    %% Condition 1: |a_1| > |a_2| & |a_3|
    [~,idx] = max([length_a1_orig,length_a2_orig,length_a3_orig]);
    lattice_vec_a(:,1)   = lattice_vec_a(:,idx);
    lattice_vec_a(:,idx) = lattice_vec_a_orig(:,1);
    Permutation([1,idx]) = Permutation([idx,1]);
    
    length_a1 = norm(lattice_vec_a(:,1));
    length_a2 = norm(lattice_vec_a(:,2));
    length_a3 = norm(lattice_vec_a(:,3));
    
    theta_1   = acos(lattice_vec_a(:,2)'*lattice_vec_a(:,3)/( length_a2*length_a3));
    theta_2   = acos(lattice_vec_a(:,1)'*lattice_vec_a(:,3)/( length_a1*length_a3));
    theta_3   = acos(lattice_vec_a(:,1)'*lattice_vec_a(:,2)/( length_a1*length_a2));
    
    %% Condition 2: ?a??sin£c_{£^}>((?a??)/(sin£c_{£^}))|cos£c_{£\}-cos£c_{£]}cos£c_{£^}|
    if length_a2*sin(theta_3) <= (length_a3/sin(theta_3))*abs(cos(theta_1)-cos(theta_2)*cos(theta_3))
        lattice_vec_a(:,[1,2,3]) = lattice_vec_a(:,[1,3,2]);
        Permutation(:,[1,2,3])   = Permutation(:,[1,3,2]);
        
        length_a1 = norm(lattice_vec_a(:,1));
        length_a2 = norm(lattice_vec_a(:,2));
        length_a3 = norm(lattice_vec_a(:,3));

        theta_1   = acos(lattice_vec_a(:,2)'*lattice_vec_a(:,3)/( length_a2*length_a3));
        theta_2   = acos(lattice_vec_a(:,1)'*lattice_vec_a(:,3)/( length_a1*length_a3));
        theta_3   = acos(lattice_vec_a(:,1)'*lattice_vec_a(:,2)/( length_a1*length_a2));
    end
    %% Test condition:  £c_£\,£c_£],£c_£^
    switch lattice_type
        case 'base_centered_monoclinic'
            if theta_1 >= pi/2
                theta_1 = pi - theta_1;
            end
            if theta_2 >= pi/2
                theta_2 = pi - theta_2;
            end
            if theta_3 >= pi/2
                theta_3 = pi - theta_3;
            end
    end
    %% Final check
    if (length_a2 > length_a1) || (length_a3 > length_a1) || length_a2*sin(theta_3) <= (length_a3/sin(theta_3))*abs(cos(theta_1)-cos(theta_2)*cos(theta_3))
        error('The lattice constants does not suitable for computation!')
    end
    
    lattice_vec_a1 = [ length_a1; 0; 0 ];
    lattice_vec_a2 = [ length_a2*cos(theta_3); length_a2*sin(theta_3); 0];
    lattice_vec_a3 = [ length_a3*cos(theta_2); 
                       length_a3*( cos(theta_1) - cos(theta_3)*cos(theta_2) )/sin(theta_3);
                       length_a3*sqrt( 1 - cos(theta_3)^2 - cos(theta_2)^2 - cos(theta_1)^2 + ...
                                  2*cos(theta_3)*cos(theta_2)*cos(theta_1) )/sin(theta_3)];
    lattice_vec_a = [lattice_vec_a1, lattice_vec_a2, lattice_vec_a3];

end

function [a1,a2,a3] = Lattice_vector_generate(lattice_type, lattice_constant)
    %% Set lattice constants
    switch lattice_type
        %% Cubic system
        case 'simple_cubic'
            a    = lattice_constant.a;
            temp = a*eye(3);
            a1   = temp(:,1);
            a2   = temp(:,2);
            a3   = temp(:,3);

        case 'face_centered_cubic'
            a    = lattice_constant.a;
%             temp = 0.5*a*[ 1, 0, 1  ;
%                            1, 1, 0  ;
%                            0, 1, 1 ];
            temp = 0.5*a*[ 0, 1, 1  ;
                           1, 0, 1  ;
                           1, 1, 0 ];
            a1   = temp(:,1);
            a2   = temp(:,2);
            a3   = temp(:,3);

        case 'body_centered_cubic'    
             a    = lattice_constant.a;
             temp = 0.5*a*[ -1,  1,  1  ;
                             1, -1,  1  ;
                             1,  1, -1 ];
             a1   = temp(:,1);
             a2   = temp(:,2);
             a3   = temp(:,3);
        %% Tetragonal system
        case 'primitive_tetragonal'
             a    = lattice_constant.a;
             c    = lattice_constant.c;
             temp = [ a, 0, 0  ;
                      0, a, 0  ;
                      0, 0, c ];
             a1   = temp(:,1);
             a2   = temp(:,2);
             a3   = temp(:,3);

        case 'body_centered_tetragonal'
             a    = lattice_constant.a;
             c    = lattice_constant.c;
             temp = 0.5*[  -a,   a,   a  ;
                            a,  -a,   a  ;
                            c,   c,  -c ];
             a1   = temp(:,1);
             a2   = temp(:,2);
             a3   = temp(:,3);   
        %% Orthorhombic system
        case 'primitive_orthorhombic'
             a    = lattice_constant.a;
             b    = lattice_constant.b;
             c    = lattice_constant.c;
             temp = [ a, 0, 0  ;
                      0, b, 0  ;
                      0, 0, c ];
             a1   = temp(:,1);
             a2   = temp(:,2);
             a3   = temp(:,3);

        case 'face_centered_orthorhombic'
             a    = lattice_constant.a;
             b    = lattice_constant.b;
             c    = lattice_constant.c;
             temp = 0.5*[ 0, a, a  ;
                          b, 0, b  ;
                          c, c, 0 ];
             a1   = temp(:,1);
             a2   = temp(:,2);
             a3   = temp(:,3);
        case 'body_centered_orthorhombic'
             a    = lattice_constant.a;
             b    = lattice_constant.b;
             c    = lattice_constant.c;
             temp = 0.5*[ -a,  a,  a  ;
                           b, -b,  b  ;
                           c,  c, -c ];
             a1   = temp(:,1);
             a2   = temp(:,2);
             a3   = temp(:,3);
        case 'a_base_centered_orthorhombic'
             a    = lattice_constant.a;
             b    = lattice_constant.b;
             c    = lattice_constant.c;
             temp = [ a,      0,      0;
                      0,  0.5*b,  0.5*b;
                      0, -0.5*c,  0.5*c];
             a1   = temp(:,1);
             a2   = temp(:,2);
             a3   = temp(:,3);

        case 'c_base_centered_orthorhombic'
             a    = lattice_constant.a;
             b    = lattice_constant.b;
             c    = lattice_constant.c;
             temp = [ 0.5*a,  0.5*a, 0;
                     -0.5*b,  0.5*b, 0;
                          0,      0, c];
             a1   = temp(:,1);
             a2   = temp(:,2);
             a3   = temp(:,3);         
        %% Hexagonal system
        case 'hexagonal' 
             a    = lattice_constant.a;
             c    = lattice_constant.c;
             %gamma = lattice_constant.gamma;
             temp = [ a,        -0.5*a, 0 ;
                      0, sqrt(3)*0.5*a, 0 ;
                      0,             0, c];
             a1   = temp(:,1);
             a2   = temp(:,2);
             a3   = temp(:,3);
        %% Rhombohedral system
        case 'rhombohedral'
             a    = lattice_constant.a;
             c    = lattice_constant.c;
             temp = [            a/2,           0,         -0.5*a   ;
                      -a/(2*sqrt(3)),   a/sqrt(3), -a/(2*sqrt(3))  ;
                                 c/3,         c/3,            c/3 ]; 
             a1   = temp(:,1);
             a2   = temp(:,2);
             a3   = temp(:,3);
        %% Monoclinic system
        case 'primitive_monoclinic'
             a    = lattice_constant.a;
             b    = lattice_constant.b;
             c    = lattice_constant.c;
             alpha = lattice_constant.alpha;
             temp = [ a,  0,             0  ;
                      0,  b,  c*cos(alpha)  ;
                      0,  0,  c*sin(alpha) ];
             a1   = temp(:,1);
             a2   = temp(:,2);
             a3   = temp(:,3);

        case 'base_centered_monoclinic'
             a     = lattice_constant.a;
             b     = lattice_constant.b;
             c     = lattice_constant.c;
             alpha = lattice_constant.alpha;
             temp  = [ 0.5*a , -0.5*a ,            0  ;
                       0.5*b ,  0.5*b , c*cos(alpha)  ;
                           0 ,      0 , c*sin(alpha) ]; 
             a1 = temp(:,1);
             a2 = temp(:,2);
             a3 = temp(:,3);
        %% Triclinic system
        case 'triclinic'    
            a     = lattice_constant.a;
            b     = lattice_constant.b;
            c     = lattice_constant.c;
            alpha = lattice_constant.alpha;
            beta  = lattice_constant.beta;
            gamma = lattice_constant.gamma;
            temp  = [ a ,  b*cos(gamma) , c*cos(beta);
                      0 ,  b*sin(gamma) , c*(cos(alpha)-cos(beta)*cos(gamma)) / sin(gamma);
                      0 ,             0 , c*( sqrt(1-cos(alpha)^2-cos(beta)^2-cos(gamma)^2+2*cos(alpha)*cos(beta)*cos(gamma)) / sin(gamma))];
            a1 = temp(:,1);
            a2 = temp(:,2);
            a3 = temp(:,3);  

    end
end