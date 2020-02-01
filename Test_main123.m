clear
clc
Nd = 100;
% A = eye(100,100);
A = rand(Nd,Nd);
A = (A' + A)/2;
        no_restart     = max(20, 2*3);
        stop_tolerance = 1.0e-6 * sqrt(Nd); %1.0e-12;
        initial_vec  = rand(Nd,1);
        initial_vec  = initial_vec / norm(initial_vec);
        target       = 0;
        target_type    = 'RGTR';%CL
        CorEq          = 'SIRA';
        LSinfo.precond = 'no';
        LSinfo.solver  = 'minres'; %'minres'; %'pcg';
        fid        = fopen('runsh_parameter.txt','r');
        [ew, ev, initial_vec, ritz_ew] = NullSp_SEP_Herm_JDSIRA_Driver_v2 (A, Nd, no_restart, 3, ...
                stop_tolerance, initial_vec, target, fid, target_type, CorEq, ...
                LSinfo,[]);%%spdiags( Lambdas.Sigma_r ,0, Nd,Nd));
        fclose(fid);
        [EV,EW] = eig(A);
        EW = diag(EW);