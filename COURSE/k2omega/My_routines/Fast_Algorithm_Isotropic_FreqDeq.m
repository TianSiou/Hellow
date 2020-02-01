function [mu, nu, Ele_field, res, Info] = Fast_Algorithm_Isotropic_FreqDeq(C, C_1, C_2, C_3, K_tilde, K_1_tilde, K_2_tilde, K_3_tilde, B, Lambdas, freq, N, Par_mesh, Par_lattice, theta, Info_LS, Info_GTSHIRA)
    % 選擇 PQEP 的 shift 值 nu_0
    Info.shift_palindromic =  exp(1i*theta); % nu_0
%     % 計算 Symplectic GEP 的 shift 值 \eta_0
%     Info.shift_ss = Info.shift_palindromic + 1/Info.shift_palindromic;
    % 調整內部線性系統的選項
    Info.linear_system.dimension   = length(C);
    Info.linear_system.tol         = 1e-12;
%     Info.linear_system.solver      = 'tfqmr';
    Info.linear_system.solver      = 'gmres';
    Info.linear_system.solver_pretreat = 'tfqmr';
    Info.linear_system.v0          = rand(3*N,1);
    Info.linear_system.isprecond   = 1;
    Info.linear_system.maxit       = 10000;
    Info.linear_system.restartSize = 30;
    % 選擇線性系統 Qp*x = b 的 precondition 參數 a_eps
    Info.linear_system.alpha = 0.25*(Info.shift_palindromic-1)^2*freq^2*sum(B.B_eps)/length(B.B_eps);
    % 調整GTSHIRA選項
    Info.GTSHIRA.matrix_size = Info.linear_system.dimension;
    Info.GTSHIRA.shift       = Info.shift_palindromic;
    Info.GTSHIRA.eigenwanted = 1;
    Info.GTSHIRA.tol         = 1e-10;
    Info.GTSHIRA.maxit       = 3;
    Info.GTSHIRA.p           = 10;
    Info.GTSHIRA.v0          = rand(2*Info.GTSHIRA.matrix_size,1);
    Info.GTSHIRA.v0          = Info.GTSHIRA.v0/norm(Info.GTSHIRA.v0);
    Info.GTSHIRA.restart_info  = 'no';
    % 比較不同 a_eps 的效果
%     NumericalResult_compare_precondition
    
    % Construst the coefficient matrices of GQEP, PQEP, and their transformations
    mtx = Matrix_Generate_1(C, C_1, C_2, C_3, K_tilde, K_1_tilde, K_2_tilde, K_3_tilde, B, Lambdas, N);
    mtx = Matrix_Generate_2(mtx, freq, Info.shift_palindromic, Info.linear_system.alpha );    
    % 建構各個 function handle
    funhand = Funhand_Generate(C, K_tilde,B, mtx, Par_mesh, Par_lattice, Lambdas, freq, Info, N);
    % 在小矩陣時畫出各種轉化下的 Spectrum
%     Test_Dense_Eigenproblems(mtx.mtx_gyroscopic, mtx.mtx_SHH, mtx.mtx_palindromic, mtx.mtx_symplecitc, mtx.mtx_ss, mtx.mtx_SH);       
    
    % 計算 GQEP & PQEP 的特徵對
    [mu, nu, Ele_field, Info] = compute_by_GTSHIRA(mtx, funhand, funhand.fun_invPp, Info);
    
    % 計算 GQEP & PQEP 特徵對的相對誤差與絕對誤差
    frobnorm_K = norm( -K_tilde.'*K_tilde, 'fro' );
    frobnorm_G = norm( C.'*K_tilde - K_tilde.'*C, 'fro' );
    frobnorm_M = norm( C.'*C - (freq^2)*spdiags(B.B_eps,0,3*N,3*N), 'fro' );
    frobnorm_A  = norm( 0.25*(C.'*C - (freq^2)*spdiags(B.B_eps,0,3*N,3*N) - (C.'*K_tilde - K_tilde.'*C) + (-K_tilde.'*K_tilde)), 'fro' );
    frobnorm_AT = norm( 0.25*(C.'*C - (freq^2)*spdiags(B.B_eps,0,3*N,3*N) + (C.'*K_tilde - K_tilde.'*C) + (-K_tilde.'*K_tilde)), 'fro' );
    frobnorm_Q  = norm( 0.5 *(C.'*C - (freq^2)*spdiags(B.B_eps,0,3*N,3*N) - (-K_tilde.'*K_tilde)), 'fro' );

    temp_g = @(lambda) abs(lambda)^2*frobnorm_M  + abs(lambda)*frobnorm_G + frobnorm_K;
    temp_p = @(lambda) abs(lambda)^2*frobnorm_AT - abs(lambda)*frobnorm_Q + frobnorm_A;
    res = {};
    for i = 1:length(mu)
        res.absres_p(i) = norm(funhand.palindromic.Qp(Ele_field(:,i),nu(i)));
        res.relres_p(i) = res.absres_p(i)/(temp_p(nu(i))*norm(Ele_field(:,i)));
        res.absres_g(i) = norm(funhand.palindromic.Qg(Ele_field(:,i),mu(i)));
        res.relres_g(i) = res.absres_g(i)/(temp_g(mu(i))*norm(Ele_field(:,i)));
    end
end

function [mu, nu, Ele_field, Info] = compute_by_GTSHIRA(mtx, funhand, fun_invPp, option)
% 用 GTSHIRA 計算 PQEP 的特徵對
    global inner_it inner_relres inner_time inner_it_pre inner_relres_pre inner_time_pre
    inner_it     = []; % initialize inner iteration number
    inner_relres = [];
    inner_time   = [];
    inner_it_pre = [];
    inner_relres_pre = [];
    inner_time_pre = [];
    Info = option;
    
    solve_LS_PQEP = @(vec_x, transp_flag) Matrix_Vector_Production_invQp( vec_x, transp_flag, mtx.mtx_palindromic, funhand, fun_invPp, option);
 
    ts = tic;
    [ nu, ev_new, outer_it ] = TSymplecticPair_GTSHIRA( option.GTSHIRA.matrix_size, funhand.palindromic.A, funhand.palindromic.Q, solve_LS_PQEP, option.GTSHIRA.shift, option.GTSHIRA.eigenwanted, option.GTSHIRA.restart_info, option.GTSHIRA );
    Info.GTSHIRA.excute_time = toc(ts);
    
    mu  = (nu - 1)./(nu + 1);
    Ele_field = ev_new(1:end/2,:);
    
    
    
    Info.GTSHIRA.outer_it     = outer_it;
    Info.linear_system.it     = inner_it;
    Info.linear_system.relres = inner_relres;
    Info.linear_system.time   = inner_time;
    Info.linear_system.it_pre     = inner_it_pre;
    Info.linear_system.relres_pre = inner_relres_pre;
    Info.linear_system.time_pre   = inner_time_pre;
end

% function [mu, nu, Ele_field ] = compute_by_eigs_on_invN_hat_K_hat( fun_invN_hat_K_hat, N, shift_palindromic )
%     % 用 eigs 計算 invN_hat_K_hat SEP 的特徵對
%     opt_eigs.isreal = 1;
%     global iter_fun_invN_hat
%     iter_fun_invN_hat = 0;
%     [z, Eta_hat, flag] = eigs(fun_invN_hat_K_hat, 6*N,2,'lm',opt_eigs);
% 
%     % 還原特徵值
%     n_ew  = length(Eta_hat);
%     Eta_0 = shift_palindromic + (1/shift_palindromic);
%     Eta   = (1./diag(Eta_hat)) + Eta_0; 
%     nu  = [0.5*(Eta + sqrt(Eta.^2-4));0.5*(Eta - sqrt(Eta.^2-4))];
%     mu  = (nu - 1)./(nu + 1);
% 
%     % 還原特徵向量
%     y_S1  = z(end/2+1:end,:); y_S2 = -z(1:end/2,:);
%     for i = 1:n_ew
%         Ele_field(:,i     ) = y_S1(:,i) + nu(i+n_ew)*y_S2(:,i);
%         Ele_field(:,i+n_ew) = y_S1(:,i) + nu(i)     *y_S2(:,i);
%     end
% end
% function [mu, nu, Ele_field ] = compute_by_eigs_on_W_hat( fun_W_hat, fun_invN2_hat, N, shift_palindromic )
%     % 用 eigs 計算 invN_hat_K_hat SEP 的特徵對
%     opt_eigs.isreal = 1;
%     global iter_fun_W_hat
%     iter_fun_W_hat = 0;
%     [z_hat, Eta_hat, flag] = eigs(fun_W_hat, 6*N,6,'lm',opt_eigs);
%     
%     for i = 1:size(z_hat,2)
%         z(:,i) = fun_invN2_hat(z_hat(:,i));
%     end
%     
%     % 還原特徵值
%     n_ew  = length(Eta_hat);
%     Eta_0 = shift_palindromic + (1/shift_palindromic);
%     Eta   = (1./diag(Eta_hat)) + Eta_0; 
%     nu  = [0.5*(Eta + sqrt(Eta.^2-4));0.5*(Eta - sqrt(Eta.^2-4))];
%     mu  = (nu - 1)./(nu + 1);
% 
%     % 還原特徵向量
%     y_S1  = z(end/2+1:end,:); y_S2 = -z(1:end/2,:);
%     for i = 1:n_ew
%         Ele_field(:,i     ) = y_S1(:,i) + nu(i+n_ew)*y_S2(:,i);
%         Ele_field(:,i+n_ew) = y_S1(:,i) + nu(i)     *y_S2(:,i);
%     end
% end
% 
% function [mu, nu, Ele_field] = compute_by_TSHIRA(mtx, N, Par_mesh, Par_lattice,freq, B, Lambdas, option)
%     % 用 TSHIRA 計算 (K_hat, N_hat) GEP 的特徵對
%     nn            = 3*N;
%     mtx_A0        = -mtx.mtx_palindromic.Q;
%     mtx_A1        =  mtx.mtx_palindromic.A.';
%     mtx_A1T       =  mtx.mtx_palindromic.A;
%     tau           = option.shift_palindromic;
%     eigswanted    = 6;
%     opts.tol      = 1e-10;
%     opts.p        = 3*eigswanted;
%     opts.maxit    = 300;
%     opts.v0       = rand(2*nn,1);
%     fun_invQp = @(sigma, vec_x, transp_flag) Matrix_Vector_Production_invQp(vec_x, transp_flag, Par_mesh.grid_num, Par_lattice.lattice_type, mtx.mtx_palindromic, Lambdas, option);
%     linear_solver_A1A0A1T = @(sigma, vec_x) fun_invQp(vec_x, 'normal');
%     linear_solver_A1TA0A1 = @(sigma, vec_x) fun_invQp(vec_x, 'trans');
%     [ev, ew] = SHIRA_TSymplP(mtx_A0, mtx_A1, mtx_A1T, linear_solver_A1A0A1T, linear_solver_A1TA0A1, nn, tau, eigswanted, opts);
%     
%     nu = ew;
%     mu  = (nu - 1)./(nu + 1);
%     Ele_field = ev;
% end