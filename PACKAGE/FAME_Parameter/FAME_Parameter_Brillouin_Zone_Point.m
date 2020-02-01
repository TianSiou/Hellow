function [ vertex ] = FAME_Parameter_Brillouin_Zone_Point(lattice_type, lattice_constant, reciprocal_lattice_vector_b)

% b = inv([a1, a2, a3]);
% reciprocal_lattice_vector_b = 2*pi*(inv(lattice_vec_a'));
% I = eye(3);
% if isfield(lattice_constant,'pretreatment_idx')
%     P = I(:,lattice_constant.pretreatment_idx);
% else
%     P = I;
% end
% reciprocal_lattice_vector_b = inv([a1, a2, a3])'*P';

switch lattice_type
    case 'simple_cubic'
        vertex.G = reciprocal_lattice_vector_b*[   0,   0,   0 ]';
        vertex.M = reciprocal_lattice_vector_b*[ 1/2, 1/2,   0 ]';
        vertex.R = reciprocal_lattice_vector_b*[ 1/2, 1/2, 1/2 ]';
        vertex.X = reciprocal_lattice_vector_b*[   0, 1/2,   0 ]';
        %============== The first Brillouin zone of a simple cubic lattice ==============
        %                            G: [   0,   0,   0 ]'   (Gamma)
        %                            M: [ 1/2, 1/2,   0 ]'
        %                            R: [ 1/2, 1/2, 1/2 ]'
        %                            X: [   0, 1/2,   0 ]'
        %================================================================================
        
    case 'face_centered_cubic'
        vertex.G = reciprocal_lattice_vector_b*[ 0  , 0  , 0   ]';
        vertex.K = reciprocal_lattice_vector_b*[ 3/8, 3/8, 3/4 ]';
        vertex.L = reciprocal_lattice_vector_b*[ 1/2, 1/2, 1/2 ]';
        vertex.U = reciprocal_lattice_vector_b*[ 5/8, 1/4, 5/8 ]';
        vertex.W = reciprocal_lattice_vector_b*[ 1/2, 1/4, 3/4 ]';
        vertex.X = reciprocal_lattice_vector_b*[ 1/2,   0, 1/2 ]';
        %=============== The first Brillouin zone of a face centered cubic lattice ==============
        %                            G = [ 0  , 0  , 0   ]'   (Gamma)
        %                            K = [ 3/8, 3/8, 3/4 ]'
        %                            L = [ 1/2, 1/2, 1/2 ]'
        %                            U = [ 5/8, 1/4, 5/8 ]'
        %                            W = [ 1/2, 1/4, 3/4 ]'
        %                            X = [ 1/2,   0, 1/2 ]'
        %=========================================================================================
        
    case 'body_centered_cubic'
        %         vertex.G  = reciprocal_lattice_vector_b*[    0,   0,   0 ]';
        %         vertex.N  = reciprocal_lattice_vector_b*[    0,   0, 1/2 ]';
        %         vertex.H  = reciprocal_lattice_vector_b*[  1/2,-1/2, 1/2 ]';
        %         vertex.P  = reciprocal_lattice_vector_b*[  1/4, 1/4, 1/4 ]';
        
        vertex.G  = reciprocal_lattice_vector_b*[    0,   0,   0 ]';
        vertex.N  = reciprocal_lattice_vector_b*[  1/2,   0,-1/2 ]';
        vertex.H  = reciprocal_lattice_vector_b*[  1/2,-1/2, 1/2 ]';
        vertex.P  = reciprocal_lattice_vector_b*[  3/4,-1/4,-1/4 ]';
        vertex.M  = reciprocal_lattice_vector_b*[  1/2, 1/2,-1/2 ]';
        %=============== The first Brillouin zone of a body centered cubic lattice ===============
        %                            G  = [    0,   0,   0 ]';   (Gamma)
        %                            N  = [    0,   0, 1/2 ]';
        %                            H  = [  1/2,-1/2, 1/2 ]';
        %                            P  = [  1/4, 1/4, 1/4 ]';
        %=========================================================================================
    case 'primitive_tetragonal'                                                              %1.4 no problem
%         vertex.G  = reciprocal_lattice_vector_b*[    0,   0,   0 ]';
%         vertex.A  = reciprocal_lattice_vector_b*[  1/2, 1/2, 1/2 ]';
%         vertex.M  = reciprocal_lattice_vector_b*[  1/2, 1/2,   0 ]';
%         vertex.R  = reciprocal_lattice_vector_b*[    0, 1/2, 1/2 ]';
%         vertex.X  = reciprocal_lattice_vector_b*[    0, 1/2,   0 ]';
%         vertex.Z  = reciprocal_lattice_vector_b*[    0,   0, 1/2 ]';
        %=============== The first Brillouin zone of a primitive_tetragonal ===============
        %                            G  = [   0,   0,   0 ]';   (Gamma)
        %                            X  = [   0, 1/2,   0 ]';
        %                            M  = [ 1/2, 1/2,   0 ]';
        %                            Z  = [   0,   0, 1/2 ]';
        %                            R  = [   0, 1/2, 1/2 ]';
        %                            A  = [ 1/2, 1/2, 1/2 ]';
        %==================================================================================
        vertex.G = reciprocal_lattice_vector_b*[   0,   0,   0 ]';
        vertex.M = reciprocal_lattice_vector_b*[ 1/2, 1/2,   0 ]';
        vertex.R = reciprocal_lattice_vector_b*[ 1/2, 1/2, 1/2 ]';
        vertex.X = reciprocal_lattice_vector_b*[   0, 1/2,   0 ]';
    case 'body_centered_tetragonal'
        if lattice_constant.c < lattice_constant.a
            eta = (1 + (lattice_constant.c/lattice_constant.a)^2)/4;
            vertex.G   = reciprocal_lattice_vector_b*[    0,     0,    0 ]';
            vertex.X   = reciprocal_lattice_vector_b*[    0,     0,  1/2 ]';
            vertex.M   = reciprocal_lattice_vector_b*[ -1/2,   1/2,  1/2 ]';
            vertex.N   = reciprocal_lattice_vector_b*[    0,   1/2,    0 ]';
            vertex.Z   = reciprocal_lattice_vector_b*[  eta,   eta, -eta ]';
            vertex.B  = reciprocal_lattice_vector_b*[ -eta, 1-eta,  eta ]';
            vertex.P   = reciprocal_lattice_vector_b*[  1/4,   1/4,  1/4 ]';
            %=============== The first Brillouin zone of a body_centered_tetragonal lattice (c < a) ===============
            % eta = (1 + (c/a)^2)/4;
            %                            G      = [     0,    0,    0 ]';   (Gamma)
            %                            X      = [    0,     0,  1/2 ]';
            %                            M      = [ -1/2,   1/2,  1/2 ]';
            %                            N      = [    0,   1/2,    0 ]';
            %                            Z      = [  eta,   eta, -eta ]';
            %                            B(Z1)  = [ -eta, 1-eta,  eta ]';
            %                            P      = [  1/4,   1/4,  1/4 ]';
            %=======================================================================================================
            
        elseif lattice_constant.c > lattice_constant.a
            eta = (1 + (lattice_constant.a/lattice_constant.c)^2)/4;
            varsigma = (lattice_constant.a^2)/(2*lattice_constant.c^2);
            vertex.G  = reciprocal_lattice_vector_b*[         0,        0,         0 ]';
            vertex.X  = reciprocal_lattice_vector_b*[         0,        0,       1/2 ]';
            vertex.Y  = reciprocal_lattice_vector_b*[ -varsigma, varsigma,       1/2 ]';
            vertex.A  = reciprocal_lattice_vector_b*[       1/2,      1/2, -varsigma ]';
            vertex.s  = reciprocal_lattice_vector_b*[      -eta,      eta,       eta ]';
            vertex.t  = reciprocal_lattice_vector_b*[       eta,    1-eta,      -eta ]';
            vertex.N  = reciprocal_lattice_vector_b*[         0,      1/2,         0 ]';
            vertex.P  = reciprocal_lattice_vector_b*[       1/4,      1/4,       1/4 ]';
            vertex.Z  = reciprocal_lattice_vector_b*[       1/2,      1/2,      -1/2 ]';
            %=============== The first Brillouin zone of a body_centered_tetragonal lattice (c > a) ===============
            % eta = (1 + (c/a)^2)/4;  varsigma = (c^2)/(2*a^2);
            %                            G     = [         0,        0,         0 ]';   (Gamma)
            %                            X     = [         0,        0,       1/2 ]';
            %                            Y     = [ -varsigma, varsigma,       1/2 ]';
            %                            A(Y1)    = [       1/2,      1/2, -varsigma ]';
            %                            s(Sig)   = [      -eta,      eta,       eta ]';   (Sigma)
            %                            t(Sig1)  = [       eta,    1-eta,      -eta ]';   (Sigma1)
            %                            N     = [         0,      1/2,         0 ]';
            %                            P     = [       1/4,      1/4,       1/4 ]';
            %                            Z     = [       1/2,      1/2,      -1/2 ]';
            %=======================================================================================================
        end
    case 'primitive_orthorhombic'
        vertex.G  = reciprocal_lattice_vector_b*[    0,   0,   0 ]';
        vertex.R  = reciprocal_lattice_vector_b*[  1/2, 1/2, 1/2 ]';
        vertex.S  = reciprocal_lattice_vector_b*[  1/2, 1/2,   0 ]';
        vertex.T  = reciprocal_lattice_vector_b*[    0, 1/2, 1/2 ]';
        vertex.U  = reciprocal_lattice_vector_b*[  1/2,   0, 1/2 ]';
        vertex.X  = reciprocal_lattice_vector_b*[  1/2,   0,   0 ]';
        vertex.Y  = reciprocal_lattice_vector_b*[    0, 1/2,   0 ]';
        vertex.Z  = reciprocal_lattice_vector_b*[    0,   0, 1/2 ]';
        
        %=============== The first Brillouin zone of a primitive_orthorhombic ===============
        %                            G  = [   0,   0,   0 ]';   (Gamma)
        %                            R  = [ 1/2, 1/2, 1/2 ]';
        %                            S  = [ 1/2, 1/2,   0 ]';
        %                            T  = [   0, 1/2, 1/2 ]';
        %                            U  = [ 1/2,   0, 1/2 ]';
        %                            X  = [ 1/2,   0,   0 ]';
        %                            Y  = [   0, 1/2,   0 ]';
        %                            Z  = [   0,   0, 1/2 ]';
        %=====================================================================================
    case 'face_centered_orthorhombic'
        if 1/(lattice_constant.a^2) >= 1/(lattice_constant.b^2)+1/(lattice_constant.c^2)
            varsigma = (1 + (lattice_constant.a/lattice_constant.b)^2 - (lattice_constant.a/lattice_constant.c)^2)/4;
            eta = (1 + (lattice_constant.a/lattice_constant.b)^2 + (lattice_constant.a/lattice_constant.c)^2)/4;
            vertex.G  = reciprocal_lattice_vector_b*[   0,            0,          0 ]';
            vertex.L  = reciprocal_lattice_vector_b*[ 1/2,          1/2,        1/2 ]';
            vertex.A  = reciprocal_lattice_vector_b*[ 1/2, 1/2+varsigma,   varsigma ]';
            vertex.B  = reciprocal_lattice_vector_b*[ 1/2, 1/2-varsigma, 1-varsigma ]';
            vertex.T  = reciprocal_lattice_vector_b*[   1,          1/2,        1/2 ]';
            vertex.Y  = reciprocal_lattice_vector_b*[ 1/2,            0,        1/2 ]';
            vertex.X  = reciprocal_lattice_vector_b*[   0,          eta,        eta ]';
            vertex.J  = reciprocal_lattice_vector_b*[   1,        1-eta,      1-eta ]';
            vertex.Z  = reciprocal_lattice_vector_b*[ 1/2,          1/2,          0 ]';
            %=============== The first Brillouin zone of a face_centered_orthorhombic lattice (1/c^2 >= 1/a^2 + 1/b^2) ===============
            % varsigma = (1 + (c/b)^2 - (c/a)^2)/4;  eta = (1 + (c/b)^2 + (c/a)^2)/4;
            %                            G      = [   0,            0,          0 ]';   (Gamma)
            %                            L      = [ 1/2,          1/2,        1/2 ]';
            %                            A      = [ 1/2, 1/2+varsigma,   varsigma ]';
            %                            B(A1)  = [ 1/2, 1/2-varsigma, 1-varsigma ]';
            %                            T      = [   1,          1/2,        1/2 ]';
            %                            Y      = [ 1/2,            0,        1/2 ]';
            %                            X      = [   0,          eta,        eta ]';
            %                            J(X1)  = [   1,        1-eta,      1-eta ]';
            %                            Z      = [ 1/2,          1/2,          0 ]';
            %========================================================================================================================
            
        else
            eta = (1 + (lattice_constant.a/lattice_constant.b)^2 - (lattice_constant.a/lattice_constant.c)^2)/4;
            delta = (1 + (lattice_constant.b/lattice_constant.a)^2 - (lattice_constant.b/lattice_constant.c)^2)/4;
            phi = (1 + (lattice_constant.c/lattice_constant.b)^2 - (lattice_constant.c/lattice_constant.a)^2)/4;
            vertex.G  = reciprocal_lattice_vector_b*[         0,       0,       0 ]';
            vertex.L  = reciprocal_lattice_vector_b*[       1/2,     1/2,     1/2 ]';
            vertex.C  = reciprocal_lattice_vector_b*[       1/2, 1/2-eta,   1-eta ]';
            vertex.E  = reciprocal_lattice_vector_b*[       1/2, 1/2+eta,     eta ]';
            vertex.D  = reciprocal_lattice_vector_b*[ 1/2-delta,     1/2, 1-delta ]';
            vertex.F  = reciprocal_lattice_vector_b*[ 1/2+delta,     1/2,   delta ]';
            vertex.H  = reciprocal_lattice_vector_b*[     1-phi, 1/2-phi,     1/2 ]';
            vertex.I  = reciprocal_lattice_vector_b*[       phi, 1/2+phi,     1/2 ]';
            vertex.X  = reciprocal_lattice_vector_b*[         0,     1/2,     1/2 ]';
            vertex.Y  = reciprocal_lattice_vector_b*[       1/2,       0,     1/2 ]';
            vertex.Z  = reciprocal_lattice_vector_b*[       1/2,     1/2,       0 ]';
            %=============== The first Brillouin zone of a face_centered_orthorhombic lattice (1/c^2 < 1/a^2 + 1/b^2) ===============
            % eta = (1 + (c/b)^2 - (c/a)^2)/4;  delta = (1 + (b/c)^2 - (b/a)^2)/4;  phi = (1 + (a/b)^2 - (a/c)^2)/4;
            %                            G      = [         0,       0,       0 ]';   (Gamma)
            %                            L      = [       1/2,     1/2,     1/2 ]';
            %                            C      = [       1/2, 1/2-eta,   1-eta ]';
            %                            E(C1)  = [       1/2, 1/2+eta,     eta ]';
            %                            D      = [ 1/2-delta,     1/2, 1-delta ]';
            %                            F(D1)  = [ 1/2+delta,     1/2,   delta ]';
            %                            H      = [     1-phi, 1/2-phi,     1/2 ]';
            %                            I(H1)  = [       phi, 1/2+phi,     1/2 ]';
            %                            X      = [         0,     1/2,     1/2 ]';
            %                            Y      = [       1/2,       0,     1/2 ]';
            %                            Z      = [       1/2,     1/2,       0 ]';
            %========================================================================================================================
        end
        
    case 'body_centered_orthorhombic'
        varsigma = (1 + (lattice_constant.a/lattice_constant.c)^2)/4;
        delta = (lattice_constant.b^2 - lattice_constant.a^2)/(4*lattice_constant.c^2);
        eta = (1 + (lattice_constant.b/lattice_constant.c)^2)/4;
        mu = (lattice_constant.a^2 +lattice_constant.b^2)/(4*lattice_constant.c^2);
        vertex.G  = reciprocal_lattice_vector_b*[         0,          0,         0 ]';
        vertex.L  = reciprocal_lattice_vector_b*[       -mu,         mu, 1/2-delta ]';
        vertex.I  = reciprocal_lattice_vector_b*[        mu,        -mu, 1/2+delta ]';
        vertex.J  = reciprocal_lattice_vector_b*[ 1/2-delta,  1/2+delta,       -mu ]';
        vertex.R  = reciprocal_lattice_vector_b*[         0,        1/2,         0 ]';
        vertex.S  = reciprocal_lattice_vector_b*[       1/2,          0,         0 ]';
        vertex.T  = reciprocal_lattice_vector_b*[         0,          0,       1/2 ]';
        vertex.W  = reciprocal_lattice_vector_b*[       1/4,        1/4,       1/4 ]';
        vertex.X  = reciprocal_lattice_vector_b*[ -varsigma,   varsigma,  varsigma ]';
        vertex.A  = reciprocal_lattice_vector_b*[  varsigma, 1-varsigma, -varsigma ]';
        vertex.Y  = reciprocal_lattice_vector_b*[       eta,       -eta,       eta ]';
        vertex.B  = reciprocal_lattice_vector_b*[     1-eta,        eta,      -eta ]';
        vertex.Z  = reciprocal_lattice_vector_b*[       1/2,        1/2,      -1/2 ]';
        %=============== The first Brillouin zone of a body_centered_orthorhombic ===============
        % varsigma = (1 + (c/a)^2)/4;  delta = (b^2 - c^2)/(4*a^2);  eta = (1 + (b/a)^2)/4;  mu = (c^2 +b^2)/(4*a^2);
        %                            G      = [         0,          0,         0 ]';   (Gamma)
        %                            L      = [       -mu,         mu, 1/2-delta ]';
        %                            I(L1)  = [        mu,        -mu, 1/2+delta ]';
        %                            J(L2)  = [ 1/2-delta,  1/2+delta,       -mu ]';
        %                            R      = [         0,        1/2,         0 ]';
        %                            S      = [       1/2,          0,         0 ]';
        %                            T      = [         0,          0,       1/2 ]';
        %                            W      = [       1/4,        1/4,       1/4 ]';
        %                            X      = [ -varsigma,   varsigma,  varsigma ]';
        %                            A(X1)  = [  varsigma, 1-varsigma, -varsigma ]';
        %                            Y      = [       eta,       -eta,       eta ]';
        %                            B(Y1)  = [     1-eta,        eta,      -eta ]';
        %                            Z      = [       1/2,        1/2,      -1/2 ]';
        %==========================================================================================
    case 'c_base_centered_orthorhombic'
        varsigma = (1 + (lattice_constant.a/lattice_constant.b)^2)/4;
        vertex.G  = reciprocal_lattice_vector_b*[         0,          0,   0 ]';
        vertex.R  = reciprocal_lattice_vector_b*[         0,        1/2, 1/2 ]';
        vertex.S  = reciprocal_lattice_vector_b*[         0,        1/2,   0 ]';
        vertex.T  = reciprocal_lattice_vector_b*[      -1/2,        1/2, 1/2 ]';
        vertex.A  = reciprocal_lattice_vector_b*[  varsigma,   varsigma, 1/2 ]';
        vertex.B  = reciprocal_lattice_vector_b*[ -varsigma, 1-varsigma, 1/2 ]';
        vertex.X  = reciprocal_lattice_vector_b*[  varsigma,   varsigma,   0 ]';
        vertex.C  = reciprocal_lattice_vector_b*[ -varsigma, 1-varsigma,   0 ]';
        vertex.Y  = reciprocal_lattice_vector_b*[      -1/2,        1/2,   0 ]';
        vertex.Z  = reciprocal_lattice_vector_b*[         0,          0, 1/2 ]';
        %=============== The first Brillouin zone of a c_base_centered_orthorhombic ===============
        % varsigma = (1 + (c/b)^2)/4;
        %                            G      = [         0,        0,     0 ]';   (Gamma)
        %                            R      = [         0,        1/2, 1/2 ]';
        %                            S      = [         0,        1/2,   0 ]';
        %                            T      = [      -1/2,        1/2, 1/2 ]';
        %                            A      = [  varsigma,   varsigma, 1/2 ]';
        %                            B(A1)  = [ -varsigma, 1-varsigma, 1/2 ]';
        %                            X      = [  varsigma,   varsigma,   0 ]';
        %                            C(X1)  = [ -varsigma, 1-varsigma,   0 ]';
        %                            Y      = [      -1/2,        1/2,   0 ]';
        %                            Z      = [         0,          0, 1/2 ]';
        %===========================================================================================
        
        
    case 'hexagonal'
        vertex.G  = reciprocal_lattice_vector_b*[    0,   0,   0 ]';
        vertex.A  = reciprocal_lattice_vector_b*[    0,   0, 1/2 ]';
        vertex.H  = reciprocal_lattice_vector_b*[  1/3, 1/3, 1/2 ]';
        vertex.K  = reciprocal_lattice_vector_b*[  1/3, 1/3,   0 ]';
        vertex.L  = reciprocal_lattice_vector_b*[  1/2,   0, 1/2 ]';
        vertex.M  = reciprocal_lattice_vector_b*[  1/2,   0,   0 ]';
        vertex.k  = reciprocal_lattice_vector_b*[  (.5*sqrt(2)-1/3), -1/3,   0 ]';
        %=============== The first Brillouin zone of a hexagonal lattice ==============
        %                            G  = [    0,   0,   0 ];   (Gamma)
        %                            A  = [    0,   0, 1/2 ];
        %                            H  = [  1/3, 1/3, 1/2 ];
        %                            K  = [  1/3, 1/3,   0 ];
        %                            L  = [  1/2,   0, 1/2 ];
        %                            M  = [  1/2,   0,   0 ];
        %==============================================================================
        
    case 'rhombohedral'
        if lattice_constant.theta_1 < pi/2
            eta = (1 + 4*cos(lattice_constant.theta_1))/(2 + 4*cos(lattice_constant.theta_1));
            nu = 3/4 - eta/2;
            vertex.G  = reciprocal_lattice_vector_b*[    0,     0,     0 ]';
            vertex.F  = reciprocal_lattice_vector_b*[  1/2,   1/2,     0 ]';
            vertex.B  = reciprocal_lattice_vector_b*[  eta,   1/2, 1-eta ]';
            vertex.C  = reciprocal_lattice_vector_b*[  1/2, 1-eta, eta-1 ]';
            vertex.L  = reciprocal_lattice_vector_b*[  1/2,     0,     0 ]';
            vertex.I  = reciprocal_lattice_vector_b*[    0,     0,  -1/2 ]';
            vertex.P  = reciprocal_lattice_vector_b*[  eta,    nu,    nu ]';
            vertex.R  = reciprocal_lattice_vector_b*[ 1-nu,  1-nu, 1-eta ]';
            vertex.S  = reciprocal_lattice_vector_b*[   nu,    nu, eta-1 ]';
            vertex.Q  = reciprocal_lattice_vector_b*[ 1-nu,    nu,     0 ]';
            vertex.X  = reciprocal_lattice_vector_b*[   nu,     0,   -nu ]';
            vertex.Z  = reciprocal_lattice_vector_b*[  1/2,   1/2,   1/2 ]';
            %=============== The first Brillouin zone of a rhombohedral lattice (theta < pi/2) ===============
            % eta = (1 + 4*cos(theta))/(2 + 4*cos(theta));  nu = 3/4 - eta/2;
            %                            G      = [    0,     0,     0 ];   (Gamma)
            %                            F      = [  1/2,   1/2,     0 ];
            %                            B      = [  eta,   1/2, 1-eta ];
            %                            C(B1)  = [  1/2, 1-eta, eta-1 ];
            %                            L      = [  1/2,     0,     0 ];
            %                            I(L1)  = [    0,     0,  -1/2 ];
            %                            P      = [  eta,    nu,    nu ];
            %                            R(P1)  = [ 1-nu,  1-nu, 1-eta ];
            %                            S(P2)  = [   nu,    nu, eta-1 ];
            %                            Q      = [ 1-nu,    nu,     0 ];
            %                            X      = [   nu,     0,   -nu ];
            %                            Z      = [  1/2,   1/2,   1/2 ];
            %=================================================================================================
            
        elseif lattice_constant.theta_1 > pi/2
            eta = 1/(2*(tan(lattice_constant.theta_1/2))^2);
            nu = 3/4 - eta/2;
            vertex.G  = reciprocal_lattice_vector_b*[     0,     0,     0 ]';
            vertex.F  = reciprocal_lattice_vector_b*[   1/2,  -1/2,     0 ]';
            vertex.L  = reciprocal_lattice_vector_b*[   1/2,     0,     0 ]';
            vertex.Z  = reciprocal_lattice_vector_b*[   1/2,  -1/2,   1/2 ]';
            vertex.P  = reciprocal_lattice_vector_b*[  1-nu,   -nu,  1-nu ]';
            vertex.R  = reciprocal_lattice_vector_b*[    nu,  nu-1,  nu-1 ]';
            vertex.Q  = reciprocal_lattice_vector_b*[   eta,   eta,   eta ]';
            vertex.S  = reciprocal_lattice_vector_b*[ 1-eta,  -eta,  -eta ]';
            %=============== The first Brillouin zone of a rhombohedral lattice (theta > pi/2) ===============
            % eta = 1/(2*tan(theta/2)^2);  nu = 3/4 - eta/2;
            %                            G      = [     0,    0,    0 ]';   (Gamma)
            %                            F      = [   1/2, -1/2,    0 ]';
            %                            L      = [   1/2,    0,    0 ]';
            %                            Z      = [   1/2, -1/2,  1/2 ]';
            %                            P      = [  1-nu,  -nu, 1-nu ]';
            %                            R(P1)  = [    nu, nu-1, nu-1 ]';
            %                            Q      = [   eta,  eta,  eta ]';
            %                            S(Q1)  = [ 1-eta, -eta,  -eta]';
            %=================================================================================================
        end
        
    case 'primitive_monoclinic'
        eta = (1 - lattice_constant.b*(cos(lattice_constant.alpha)/lattice_constant.c)) /(2*(sin(lattice_constant.alpha))^2);
        nu = 1/2 - eta*lattice_constant.c*(cos(lattice_constant.alpha)/lattice_constant.b);
        vertex.G  = reciprocal_lattice_vector_b*[   0,     0,    0 ]';
        vertex.A  = reciprocal_lattice_vector_b*[ 1/2,   1/2,    0 ]';
        vertex.C  = reciprocal_lattice_vector_b*[   0,   1/2,  1/2 ]';
        vertex.D  = reciprocal_lattice_vector_b*[ 1/2,     0,  1/2 ]';
        vertex.F  = reciprocal_lattice_vector_b*[ 1/2,     0, -1/2 ]';
        vertex.E  = reciprocal_lattice_vector_b*[ 1/2,   1/2,  1/2 ]';
        vertex.H  = reciprocal_lattice_vector_b*[   0,   eta, 1-nu ]';
        vertex.I  = reciprocal_lattice_vector_b*[   0, 1-eta,   nu ]';
        vertex.J  = reciprocal_lattice_vector_b*[   0,   eta,  -nu ]';
        vertex.M  = reciprocal_lattice_vector_b*[ 1/2,   eta, 1-nu ]';
        vertex.N  = reciprocal_lattice_vector_b*[ 1/2, 1-eta,   nu ]';
        vertex.O  = reciprocal_lattice_vector_b*[ 1/2,   eta,  -nu ]';
        vertex.X  = reciprocal_lattice_vector_b*[   0,   1/2,    0 ]';
        vertex.Y  = reciprocal_lattice_vector_b*[   0,     0,  1/2 ]';
        vertex.B  = reciprocal_lattice_vector_b*[   0,     0, -1/2 ]';
        vertex.Z  = reciprocal_lattice_vector_b*[ 1/2,     0,    0 ]';
        %=============== The first Brillouin zone of a primitive_monoclinic ===============
        % eta = (1 - b*(cos(gamma/a)))/(2*(sin(gamma))^2);  nu = 1/2 - eta*a*(cos(gamma/b));
        %                            G      = [   0,     0,    0 ]';   (Gamma)
        %                            A      = [ 1/2,   1/2,    0 ]';
        %                            C      = [   0,   1/2,  1/2 ]';
        %                            D      = [ 1/2,     0,  1/2 ]';
        %                            F(D1)  = [ 1/2,     0, -1/2 ]';
        %                            E      = [ 1/2,   1/2,  1/2 ]';
        %                            H      = [   0,   eta, 1-nu ]';
        %                            I(H1)  = [   0, 1-eta,   nu ]';
        %                            J(H2)  = [   0,   eta,  -nu ]';
        %                            M      = [ 1/2,   eta, 1-nu ]';
        %                            N(M1)  = [ 1/2, 1-eta,   nu ]';
        %                            O(M2)  = [ 1/2,   eta,  -nu ]';
        %                            X      = [   0,   1/2,    0 ]';
        %                            Y      = [   0,     0,  1/2 ]';
        %                            B(Y1)  = [   0,     0, -1/2 ]';
        %                            Z      = [ 1/2,     0,    0 ]';
        %===================================================================================
        
    case 'base_centered_monoclinic'
        if reciprocal_lattice_vector_b(:,1).'*reciprocal_lattice_vector_b(:,2) <= 0   % k_gamma >= pi/2
            varsigma = (2 - lattice_constant.b*cos(lattice_constant.alpha)/lattice_constant.c)/(4*sin(lattice_constant.alpha)^2);
            eta      = 1/2 + 2*varsigma*lattice_constant.c*cos(lattice_constant.alpha)/lattice_constant.b;
            psi      = 3/4 - (lattice_constant.a^2)/(4*lattice_constant.b^2*sin(lattice_constant.alpha)^2);
            phi      = psi + (3/4 - psi)*lattice_constant.b*cos(lattice_constant.alpha)/lattice_constant.c;
            vertex.G  = reciprocal_lattice_vector_b*[          0,          0,     0 ]';
            vertex.N  = reciprocal_lattice_vector_b*[        1/2,          0,     0 ]';
            vertex.S  = reciprocal_lattice_vector_b*[          0,       -1/2,     0 ]';
            vertex.F  = reciprocal_lattice_vector_b*[ 1-varsigma, 1-varsigma, 1-eta ]';
            vertex.J  = reciprocal_lattice_vector_b*[   varsigma,   varsigma,   eta ]';
            vertex.K  = reciprocal_lattice_vector_b*[  -varsigma,  -varsigma, 1-eta ]';
            vertex.O  = reciprocal_lattice_vector_b*[ 1-varsigma,  -varsigma, 1-eta ]';
            vertex.I  = reciprocal_lattice_vector_b*[        phi,      1-phi,   1/2 ]';
            vertex.R  = reciprocal_lattice_vector_b*[      1-phi,      phi-1,   1/2 ]';
            vertex.L  = reciprocal_lattice_vector_b*[        1/2,        1/2,   1/2 ]';
            vertex.M  = reciprocal_lattice_vector_b*[        1/2,          0,   1/2 ]';
            vertex.X  = reciprocal_lattice_vector_b*[      1-psi,      psi-1,     0 ]';
            vertex.T  = reciprocal_lattice_vector_b*[        psi,      1-psi,     0 ]';
            vertex.U  = reciprocal_lattice_vector_b*[      psi-1,       -psi,     0 ]';
            vertex.Y  = reciprocal_lattice_vector_b*[        1/2,        1/2,     0 ]';
            vertex.V  = reciprocal_lattice_vector_b*[       -1/2,       -1/2,     0 ]';
            vertex.Z  = reciprocal_lattice_vector_b*[          0,          0,   1/2 ]';
            %=============== The first Brillouin zone of a a_base_centered_monoclinic lattice (k_alpha >= pi/2) ===============
            % varsigma = (2 - b*cos(gamma)/a)/(4*sin(gamma)^2);  eta = 1/2 + 2*a*cos(gamma)/b;  psi = 3/4 - (c^2)/(4*b^2*sin(gamma)^2);
            % phi = psi + (3/4 - psi)*b*cos(gamma)/a;
            %                            G      = [          0,          0,     0 ]';   (Gamma)
            %                            N      = [        1/2,          0,     0 ]';
            %                            S(N1)  = [          0,       -1/2,     0 ]';
            %                            F      = [ 1-varsigma, 1-varsigma, 1-eta ]';
            %                            J(F1)  = [   varsigma,   varsigma,   eta ]';
            %                            K(F2)  = [  -varsigma,  -varsigma, 1-eta ]';
            %                            O(F3)  = [ 1-varsigma,  -varsigma, 1-eta ]';
            %                            I      = [        phi,      1-phi,   1/2 ]';
            %                            R(I1)  = [      1-phi,      phi-1,   1/2 ]';
            %                            L      = [        1/2,        1/2,   1/2 ]';
            %                            M      = [        1/2,          0,   1/2 ]';
            %                            X      = [      1-psi,      psi-1,     0 ]';
            %                            T(X1)  = [        psi,      1-psi,     0 ]';
            %                            U(X2)  = [      psi-1,       -psi,     0 ]';
            %                            Y      = [        1/2,        1/2,     0 ]';
            %                            V(Y1)  = [       -1/2,       -1/2,     0 ]';
            %                            Z      = [          0,          0,   1/2 ]';
            %===============================================================================================================
            
        else
            if lattice_constant.b*cos(lattice_constant.alpha)/lattice_constant.c + lattice_constant.b^2*sin(lattice_constant.alpha)^2/lattice_constant.a^2 <= 1
                mu       = (1 + (lattice_constant.b/lattice_constant.a)^2)/4;
                delta    = lattice_constant.b*lattice_constant.c*cos(lattice_constant.alpha)/(2*lattice_constant.a^2);
                varsigma = mu - .25 + (1 - lattice_constant.b*cos(lattice_constant.alpha)/lattice_constant.c)/(4*sin(lattice_constant.alpha)^2);
                eta      = .5 + 2*varsigma*lattice_constant.c*cos(lattice_constant.alpha)/lattice_constant.b;
                phi      = 1 + varsigma - 2*mu;
                psi      = eta - 2*delta;
                vertex.G  = reciprocal_lattice_vector_b*[          0,         0,      0 ]';
                vertex.F  = reciprocal_lattice_vector_b*[      1-phi,     1-phi,  1-psi ]';
                vertex.J  = reciprocal_lattice_vector_b*[        phi,     phi-1,    psi ]';
                vertex.K  = reciprocal_lattice_vector_b*[      1-phi,      -phi,  1-psi ]';
                vertex.H  = reciprocal_lattice_vector_b*[   varsigma,  varsigma,    eta ]';
                vertex.P  = reciprocal_lattice_vector_b*[ 1-varsigma, -varsigma,  1-eta ]';
                vertex.Q  = reciprocal_lattice_vector_b*[  -varsigma, -varsigma,  1-eta ]';
                vertex.I  = reciprocal_lattice_vector_b*[        1/2,      -1/2,    1/2 ]';
                vertex.M  = reciprocal_lattice_vector_b*[        1/2,         0,    1/2 ]';
                vertex.N  = reciprocal_lattice_vector_b*[        1/2,         0,      0 ]';
                vertex.S  = reciprocal_lattice_vector_b*[          0,      -1/2,      0 ]';
                vertex.X  = reciprocal_lattice_vector_b*[        1/2,      -1/2,      0 ]';
                vertex.Y  = reciprocal_lattice_vector_b*[         mu,        mu,  delta ]';
                vertex.V  = reciprocal_lattice_vector_b*[       1-mu,       -mu, -delta ]';
                vertex.W  = reciprocal_lattice_vector_b*[        -mu,       -mu, -delta ]';
                vertex.A  = reciprocal_lattice_vector_b*[         mu,      mu-1,  delta ]';
                vertex.Z  = reciprocal_lattice_vector_b*[          0,         0,    1/2 ]';
                %=============== The first Brillouin zone of a a_base_centered_monoclinic lattice ( b*cos(gamma)/a + b^2*sin(gamma)^2/c^2 <= 1) ===============
                % mu = (1 + (b/c)^2)/4;  delta = b*a*cos(gamma)/(2*c^2);  varsigma = mu - 1/4 + (1 - b*cos(gamma)/a)/(4*sin(gamma)^2);
                % eta = 1/2 + 2*varsigma*a*cos(gamma)/b;  phi = 1 + varsigma - 2*mu;  psi = eta - 2*delta;
                %                            G      = [          0,         0,      0 ]';   (Gamma)
                %                            F      = [      1-phi,     1-phi,  1-psi ]';
                %                            J(F1)  = [        phi,     phi-1,    psi ]';
                %                            K(F2)  = [      1-phi,      -phi,  1-psi ]';
                %                            H      = [   varsigma,  varsigma,    eta ]';
                %                            P(H1)  = [ 1-varsigma, -varsigma,  1-eta ]';
                %                            Q(H2)  = [  -varsigma, -varsigma,  1-eta ]';
                %                            I      = [        1/2,      -1/2,    1/2 ]';
                %                            M      = [        1/2,         0,    1/2 ]';
                %                            N      = [        1/2,         0,      0 ]';
                %                            S(N1)  = [          0,      -1/2,      0 ]';
                %                            X      = [        1/2,      -1/2,      0 ]';
                %                            Y      = [         mu,        mu,  delta ]';
                %                            V(Y1)  = [       1-mu,       -mu, -delta ]';
                %                            W(Y2)  = [        -mu,       -mu, -delta ]';
                %                            A(Y3)  = [         mu,      mu-1,  delta ]';
                %                            Z      = [          0,         0,    1/2 ]';
                %============================================================================================================================================
                
            else
                varsigma = ((lattice_constant.b/lattice_constant.a)^2 + (1 - lattice_constant.b*cos(lattice_constant.alpha)/lattice_constant.c)/sin(lattice_constant.alpha)^2)/4;
                eta = .5 + 2*varsigma*lattice_constant.c*cos(lattice_constant.alpha)/lattice_constant.b;
                mu = eta/2 + lattice_constant.b^2/(4*lattice_constant.a^2) - lattice_constant.b*lattice_constant.c*cos(lattice_constant.alpha)/(2*lattice_constant.a^2);
                nu = 2*mu - varsigma;
                omega = ((4*nu -1 - lattice_constant.b^2*sin(lattice_constant.alpha)^2/lattice_constant.a^2)*lattice_constant.c)/(2*lattice_constant.b*cos(lattice_constant.alpha));
                delta = varsigma*lattice_constant.c*cos(lattice_constant.alpha)/lattice_constant.b + omega/2 - 1/4;
                rho = 1 - varsigma*(lattice_constant.a/lattice_constant.b)^2;
                vertex.G  = reciprocal_lattice_vector_b*[          0,         0,       0 ]';
                vertex.F  = reciprocal_lattice_vector_b*[         nu,        nu,   omega ]';
                vertex.J  = reciprocal_lattice_vector_b*[       1-nu,      1-nu, 1-omega ]';
                vertex.K  = reciprocal_lattice_vector_b*[         nu,      nu-1,   omega ]';
                vertex.H  = reciprocal_lattice_vector_b*[   varsigma,  varsigma,     eta ]';
                vertex.P  = reciprocal_lattice_vector_b*[ 1-varsigma, -varsigma,   1-eta ]';
                vertex.Q  = reciprocal_lattice_vector_b*[  -varsigma, -varsigma,   1-eta ]';
                vertex.I  = reciprocal_lattice_vector_b*[        rho,     1-rho,     1/2 ]';
                vertex.R  = reciprocal_lattice_vector_b*[      1-rho,     rho-1,     1/2 ]';
                vertex.L  = reciprocal_lattice_vector_b*[        1/2,       1/2,     1/2 ]';
                vertex.M  = reciprocal_lattice_vector_b*[        1/2,         0,     1/2 ]';
                vertex.N  = reciprocal_lattice_vector_b*[        1/2,         0,       0 ]';
                vertex.S  = reciprocal_lattice_vector_b*[          0,      -1/2,       0 ]';
                vertex.X  = reciprocal_lattice_vector_b*[        1/2,      -1/2,       0 ]';
                vertex.Y  = reciprocal_lattice_vector_b*[         mu,        mu,   delta ]';
                vertex.V  = reciprocal_lattice_vector_b*[       1-mu,       -mu,  -delta ]';
                vertex.W  = reciprocal_lattice_vector_b*[        -mu,       -mu,  -delta ]';
                vertex.A  = reciprocal_lattice_vector_b*[         mu,      mu-1,   delta ]';
                vertex.Z  = reciprocal_lattice_vector_b*[          0,         0,     1/2 ]';
                %=============== The first Brillouin zone of a a_base_centered_monoclinic lattice ( b*cos(gamma)/a + b^2*sin(gamma)^2/c^2 > 1) ===============
                % varsigma = ((b/c)^2 + (1 - b*cos(gamma)/a)/sin(gamma)^2)/4;  eta = 1/2 + 2*varsigma*a*cos(gamma)/b;
                % mu = eta/2 + b^2/(4*c^2) - b*a*cos(gamma)/(2*c^2);  nu = 2*mu - varsigma;
                % omega = ((4*nu -1 - b^2*sin(gamma)^2/c^2)*a)/(2*b*cos(gamma));  delta = varsigma*a*cos(gamma)/b + omega/2 - 1/4;  rho = 1 - varsigma*(c/b)^2;
                %                            G      = [          0,         0,       0 ]';   (Gamma)
                %                            F      = [         nu,        nu,   omega ]';
                %                            J(F1)  = [       1-nu,      1-nu, 1-omega ]';
                %                            K(F2)  = [         nu,      nu-1,   omega ]';
                %                            H      = [   varsigma,  varsigma,     eta ]';
                %                            P(H1)  = [ 1-varsigma, -varsigma,   1-eta ]';
                %                            Q(H2)  = [  -varsigma, -varsigma,   1-eta ]';
                %                            I      = [        rho,     1-rho,     1/2 ]';
                %                            R(I1)  = [      1-rho,     rho-1,     1/2 ]';
                %                            L      = [        1/2,       1/2,     1/2 ]';
                %                            M      = [        1/2,         0,     1/2 ]';
                %                            N      = [        1/2,         0,       0 ]';
                %                            S(N1)  = [          0,      -1/2,       0 ]';
                %                            X      = [        1/2,      -1/2,       0 ]';
                %                            Y      = [         mu,        mu,   delta ]';
                %                            V(Y1)  = [       1-mu,       -mu,  -delta ]';
                %                            W(Y2)  = [        -mu,       -mu,  -delta ]';
                %                            A(Y3)  = [         mu,      mu-1,   delta ]';
                %                            Z      = [          0,         0,     1/2 ]';
                %============================================================================================================================================
            end
        end
        
    case 'triclinic'
        %         if reciprocal_lattice_vector_b(:,2).'*reciprocal_lattice_vector_b(:,3) < 0   % k_alpha > pi/2
        %             vertex.G  = reciprocal_lattice_vector_b*[   0,   0,   0 ]';
        %             vertex.L  = reciprocal_lattice_vector_b*[ 1/2, 1/2,   0 ]';
        %             vertex.M  = reciprocal_lattice_vector_b*[   0, 1/2, 1/2 ]';
        %             vertex.N  = reciprocal_lattice_vector_b*[ 1/2,   0, 1/2 ]';
        %             vertex.R  = reciprocal_lattice_vector_b*[ 1/2, 1/2, 1/2 ]';
        %             vertex.X  = reciprocal_lattice_vector_b*[ 1/2,   0,   0 ]';
        %             vertex.Y  = reciprocal_lattice_vector_b*[   0, 1/2,   0 ]';
        %             vertex.Z  = reciprocal_lattice_vector_b*[   0,   0, 1/2 ]';
        % %=============== The first Brillouin zone of a triclinic lattice (k_alpha > pi/2) ===============
        % %                            G   = [   0,   0,   0 ]';   (Gamma)
        % %                            L   = [ 1/2, 1/2,   0 ]';
        % %                            M   = [   0, 1/2, 1/2 ]';
        % %                            N   = [ 1/2,   0, 1/2 ]';
        % %                            R   = [ 1/2, 1/2, 1/2 ]';
        % %                            X   = [ 1/2,   0,   0 ]';
        % %                            Y   = [   0, 1/2,   0 ]';
        % %                            Z   = [   0,   0, 1/2 ]';
        % %=======================================================================
        %
        %         elseif reciprocal_lattice_vector_b(:,2).'*reciprocal_lattice_vector_b(:,3) > 0   % k_alpha < pi/2
        %             vertex.G  = reciprocal_lattice_vector_b*[    0,    0,   0 ]';
        %             vertex.L  = reciprocal_lattice_vector_b*[  1/2, -1/2,   0 ]';
        %             vertex.M  = reciprocal_lattice_vector_b*[    0,    0, 1/2 ]';
        %             vertex.N  = reciprocal_lattice_vector_b*[ -1/2, -1/2, 1/2 ]';
        %             vertex.R  = reciprocal_lattice_vector_b*[    0, -1/2, 1/2 ]';
        %             vertex.X  = reciprocal_lattice_vector_b*[    0, -1/2,   0 ]';
        %             vertex.Y  = reciprocal_lattice_vector_b*[  1/2,    0,   0 ]';
        %             vertex.Z  = reciprocal_lattice_vector_b*[ -1/2,    0, 1/2 ]';
        % %=============== The first Brillouin zone of a triclinic (k_alpha < pi/2) ===============
        % %                            G   = [    0,    0,   0 ]';   (Gamma)
        % %                            L   = [  1/2, -1/2,   0 ]';
        % %                            M   = [    0,    0, 1/2 ]';
        % %                            N   = [ -1/2, -1/2, 1/2 ]';
        % %                            R   = [    0, -1/2, 1/2 ]';
        % %                            X   = [    0, -1/2,   0 ]';
        % %                            Y   = [  1/2,    0,   0 ]';
        % %                            Z   = [ -1/2,    0, 1/2 ]';
        % %=======================================================================
        %         end
        vertex.G  = reciprocal_lattice_vector_b*[    0,    0,    0 ]'; % [     0,    0,    0]
        vertex.N  = reciprocal_lattice_vector_b*[  1/2,    0, -1/2 ]'; % [  -0.5,    0,  0.5]
        vertex.H  = reciprocal_lattice_vector_b*[  1/2, -1/2,  1/2 ]'; % [     0,    1,    0]
        vertex.P  = reciprocal_lattice_vector_b*[  3/4, -1/4, -1/4 ]'; % [  -0.5,  0.5,  0.5]
        vertex.M  = reciprocal_lattice_vector_b*[  1/2,  1/2, -1/2 ]'; % [     0,    0,    1]
        %         vertex.W  = reciprocal_lattice_vector_b*[ 0.08, 0.12, 0.08 ]'; % [   0.2, 0.16,  0.2]
        %         vertex.E  = reciprocal_lattice_vector_b*[  0.1,  0.1,  0.1 ]'; % [   0.2,  0.2,  0.2]
        %         vertex.R  = reciprocal_lattice_vector_b*[  0.1, -0.1,  0.1 ]'; % [     0,  0.2,    0]
        vertex.Q  = reciprocal_lattice_vector_b*[ -1/4, -1/4,  3/4 ]'; % [   0.5,  0.5, -0.5]
        vertex.R  = reciprocal_lattice_vector_b*[  0.7, -0.2, -0.3 ]'; % [  -0.5,  0.4,  0.5]
        vertex.S  = reciprocal_lattice_vector_b*[ -0.3, -0.2,  0.7 ]'; % [   0.5,  0.4, -0.5]
%         vertex.H  = [0, 1, 0]';
end
end