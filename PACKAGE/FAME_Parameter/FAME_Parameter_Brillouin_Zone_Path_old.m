function [ path_point, path_point_num ] = FAME_Parameter_Brillouin_Zone_Path( path_string, part_num, lattice_type, lattice_constant )
%============= The first Brillouin zone of a body centered cubic lattice ======
%                            G  = [    0,   0,   0 ];      
%                            H  = [    0,   1,   0 ];
%                            Z  = [    0,   0,   1 ];
%                            h  = [    0, 3/4,   0 ];
%                            P  = [ -1/2, 1/2, 1/2 ];
%                            p  = [ -1/2, 1/4, 1/2 ];
%                            N  = [ -1/2,   0, 1/2 ]; 
%============= The first Brillouin zone of a simple cubic lattice ======
%                            G: [   0,   0,   0 ]'
%                            X: [ 0.5,   0,   0 ]'
%                            M: [ 0.5, 0.5,   0 ]'
%                            R: [ 0.5, 0.5, 0.5 ]'
%============= The first Brillouin zone of a face centered cubic lattice ======
%                            X = [ 0  , 1  , 0   ]'          
%                            U = [ 1/4, 1  , 1/4 ]'
%                            L = [ 1/2, 1/2, 1/2 ]'
%                            G = [ 0  , 0  , 0   ]'                                                                   ;
%                            W = [ 1/2, 1  , 0   ]'
%                            K = [ 3/4, 3/4, 0   ]'
%==============================================================================
    n = length(path_string);
    switch lattice_type
        case 'simple_cubic'
            vertex.G = [  0,  0,  0 ]/lattice_constant.a;
            vertex.X = [  0, .5,  0 ]/lattice_constant.a;
            vertex.M = [ .5, .5,  0 ]/lattice_constant.a;
            vertex.R = [ .5, .5, .5 ]/lattice_constant.a;
        case 'face_centered_cubic'
           Q = (1/sqrt(2))*[ 1        , 1         , 0         ;
                            -1/sqrt(3), 1/sqrt(3) , 2/sqrt(3) ;
                             2/sqrt(6),-2/sqrt(6) , 2/sqrt(6)];
            vertex.X = [ 0  , 1  , 0   ]*Q'/lattice_constant.a;             
            vertex.U = [ 1/4, 1  , 1/4 ]*Q'/lattice_constant.a;
            vertex.L = [ 1/2, 1/2, 1/2 ]*Q'/lattice_constant.a;
            vertex.G = [ 0  , 0  , 0   ]                                                                    ;
            vertex.W = [ 1/2, 1  , 0   ]*Q'/lattice_constant.a;
            vertex.K = [ 3/4, 3/4, 0   ]*Q'/lattice_constant.a;
        case 'body_centered_cubic'
            Q = (1/sqrt(3))*[      -1   ,       1   ,     1   ;
                               sqrt(2)/2, -sqrt(2)/2, sqrt(2) ;
                               sqrt(6)/2,  sqrt(6)/2,     0  ];
            vertex.G  = [    0,   0,   0 ];      
            vertex.H  = [    0,   1,   0 ]*Q'/lattice_constant.a;
            vertex.Z  = [    0,   0,   1 ]*Q'/lattice_constant.a;
            vertex.h  = [    0, 3/4,   0 ]*Q'/lattice_constant.a;
            vertex.P  = [ -1/2, 1/2, 1/2 ]*Q'/lattice_constant.a;
            vertex.p  = [ -1/2, 1/4, 1/2 ]*Q'/lattice_constant.a;
            vertex.N  = [ -1/2,   0, 1/2 ]*Q'/lattice_constant.a; 
%         case 'general'
        otherwise
           Q = (1/sqrt(2))*[ 1        , 1         , 0         ;
                            -1/sqrt(3), 1/sqrt(3) , 2/sqrt(3) ;
                             2/sqrt(6),-2/sqrt(6) , 2/sqrt(6)];
            vertex.X = [ 0  , 1  , 0   ]*Q';             
            vertex.U = [ 1/4, 1  , 1/4 ]*Q';
            vertex.L = [ 1/2, 1/2, 1/2 ]*Q';
            vertex.G = [ 0  , 0  , 0   ]                                                                    ;
            vertex.W = [ 1/2, 1  , 0   ]*Q';
            vertex.K = [ 3/4, 3/4, 0   ]*Q';
    end
    path_point(:,1) = eval( ['vertex.',path_string(1)] );
    
    for i = 1:n-1
        switch path_string(i)
            case 'G'
                path_point = [ path_point, divide_subpath(['G',path_string(i+1)], vertex.G, part_num, vertex) ];
            case 'X'
                path_point = [ path_point, divide_subpath(['X',path_string(i+1)], vertex.X, part_num, vertex) ]; 
            case 'M'
                path_point = [ path_point, divide_subpath(['M',path_string(i+1)], vertex.M, part_num, vertex) ];
            case 'R'
                path_point = [ path_point, divide_subpath(['R',path_string(i+1)], vertex.R, part_num, vertex) ];
            case 'U'
                path_point = [ path_point, divide_subpath(['U',path_string(i+1)], vertex.U, part_num, vertex) ]; 
            case 'L'
                path_point = [ path_point, divide_subpath(['L',path_string(i+1)], vertex.L, part_num, vertex) ];
            case 'W'
                path_point = [ path_point, divide_subpath(['W',path_string(i+1)], vertex.W, part_num, vertex) ];
            case 'K'
                path_point = [ path_point, divide_subpath(['K',path_string(i+1)], vertex.K, part_num, vertex) ];
            case 'H'
                path_point = [ path_point, divide_subpath(['H',path_string(i+1)], vertex.H, part_num, vertex) ]; 
            case 'Z'
                path_point = [ path_point, divide_subpath(['Z',path_string(i+1)], vertex.Z, part_num, vertex) ]; 
            case 'h'
                path_point = [ path_point, divide_subpath(['h',path_string(i+1)], vertex.h, part_num, vertex) ]; 
            case 'P'
                path_point = [ path_point, divide_subpath(['P',path_string(i+1)], vertex.P, part_num, vertex) ];
            case 'p'
                path_point = [ path_point, divide_subpath(['p',path_string(i+1)], vertex.p, part_num, vertex) ];
            case 'N'
                path_point = [ path_point, divide_subpath(['N',path_string(i+1)], vertex.N, part_num, vertex) ];                
        end
    end
    path_point_num = size(path_point,2);
end

function subpath_point = divide_subpath(subpath_string, first_point, part_num, vertex )
    
    switch subpath_string(2)
        case 'G'
            subpath_point(1,:) = linspace(first_point(1), vertex.G(1), part_num);
            subpath_point(2,:) = linspace(first_point(2), vertex.G(2), part_num);
            subpath_point(3,:) = linspace(first_point(3), vertex.G(3), part_num);
        case 'X'
            subpath_point(1,:) = linspace(first_point(1), vertex.X(1), part_num);
            subpath_point(2,:) = linspace(first_point(2), vertex.X(2), part_num);
            subpath_point(3,:) = linspace(first_point(3), vertex.X(3), part_num);
        case 'M'
            subpath_point(1,:) = linspace(first_point(1), vertex.M(1), part_num);
            subpath_point(2,:) = linspace(first_point(2), vertex.M(2), part_num);
            subpath_point(3,:) = linspace(first_point(3), vertex.M(3), part_num);
        case 'R'
            subpath_point(1,:) = linspace(first_point(1), vertex.R(1), part_num);
            subpath_point(2,:) = linspace(first_point(2), vertex.R(2), part_num);
            subpath_point(3,:) = linspace(first_point(3), vertex.R(3), part_num);
        case 'U'
            subpath_point(1,:) = linspace(first_point(1), vertex.U(1), part_num);
            subpath_point(2,:) = linspace(first_point(2), vertex.U(2), part_num);
            subpath_point(3,:) = linspace(first_point(3), vertex.U(3), part_num);
        case 'L'
            subpath_point(1,:) = linspace(first_point(1), vertex.L(1), part_num);
            subpath_point(2,:) = linspace(first_point(2), vertex.L(2), part_num);
            subpath_point(3,:) = linspace(first_point(3), vertex.L(3), part_num);
        case 'W'
            subpath_point(1,:) = linspace(first_point(1), vertex.W(1), part_num);
            subpath_point(2,:) = linspace(first_point(2), vertex.W(2), part_num);
            subpath_point(3,:) = linspace(first_point(3), vertex.W(3), part_num);
        case 'K'
            subpath_point(1,:) = linspace(first_point(1), vertex.K(1), part_num);
            subpath_point(2,:) = linspace(first_point(2), vertex.K(2), part_num);
            subpath_point(3,:) = linspace(first_point(3), vertex.K(3), part_num);   
        case 'H'
            subpath_point(1,:) = linspace(first_point(1), vertex.H(1), part_num);
            subpath_point(2,:) = linspace(first_point(2), vertex.H(2), part_num);
            subpath_point(3,:) = linspace(first_point(3), vertex.H(3), part_num);
        case 'Z'
            subpath_point(1,:) = linspace(first_point(1), vertex.Z(1), part_num);
            subpath_point(2,:) = linspace(first_point(2), vertex.Z(2), part_num);
            subpath_point(3,:) = linspace(first_point(3), vertex.Z(3), part_num);
        case 'h'
            subpath_point(1,:) = linspace(first_point(1), vertex.h(1), part_num);
            subpath_point(2,:) = linspace(first_point(2), vertex.h(2), part_num);
            subpath_point(3,:) = linspace(first_point(3), vertex.h(3), part_num);
        case 'P'
            subpath_point(1,:) = linspace(first_point(1), vertex.P(1), part_num);
            subpath_point(2,:) = linspace(first_point(2), vertex.P(2), part_num);
            subpath_point(3,:) = linspace(first_point(3), vertex.P(3), part_num);
        case 'p'
            subpath_point(1,:) = linspace(first_point(1), vertex.p(1), part_num);
            subpath_point(2,:) = linspace(first_point(2), vertex.p(2), part_num);
            subpath_point(3,:) = linspace(first_point(3), vertex.p(3), part_num);
        case 'N'
            subpath_point(1,:) = linspace(first_point(1), vertex.N(1), part_num);
            subpath_point(2,:) = linspace(first_point(2), vertex.N(2), part_num);
            subpath_point(3,:) = linspace(first_point(3), vertex.N(3), part_num);            
    end
    subpath_point = subpath_point(:,2:end);
end