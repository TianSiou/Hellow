function NumericalResult_compare_precondition(C, C_1, C_2, C_3, K_tilde, K_1_tilde, K_2_tilde, K_3_tilde, B, Lambdas, freq, N, Par_mesh, Par_lattice, FAME_option)
% Compare the efficiency between different precondition parameter alpha and
% shift value nu_0
    
%     theta = linspace(0,2*pi,10);
%     theta = theta(2:end-1);
    theta = pi/2;
    ratio = 0.5:0.1:1.5;
    for i = 1:length(ratio)
        % 選擇 PQEP 的 shift 值 nu_0
        Info.shift_palindromic =  exp(1i*theta); % nu_0
        % 計算 Symplectic GEP 的 shift 值 \eta_0
        Info.shift_ss = Info.shift_palindromic + 1/Info.shift_palindromic;
        % 調整內部線性系統的選項
        Info.linear_system.dimension   = length(C);
        Info.linear_system.tol         = 1e-12;
        Info.linear_system.solver      = 'gmres';
        Info.linear_system.solver_pretreat = 'tfqmr';
        Info.linear_system.v0          = ones(length(C),1)/sqrt(length(C));
        Info.linear_system.isprecond   = 1;
        Info.linear_system.maxit       = 10000;
        Info.linear_system.restartSize = 30;
        % 選擇線性系統 Qp*x = b 的 precondition 參數 a_eps
        a_eps = 0.25*(Info.shift_palindromic-1)^2*freq^2*sum(B.B_eps)/length(B.B_eps);
        Info.linear_system.alpha = ratio(i)*a_eps;
        
        mtx     = Matrix_Generate_1(C, C_1, C_2, C_3, K_tilde, K_1_tilde, K_2_tilde, K_3_tilde, B, Lambdas, N);
        mtx     = Matrix_Generate_2(mtx, freq, Info.shift_palindromic, Info.linear_system.alpha );    
        
        funhand = Funhand_Generate(C, K_tilde,B, mtx, Par_mesh, Par_lattice, Lambdas, freq, Info, N);

        vec_x = ones(length(C),1)/sqrt(length(C));
        [ vec_y  , relres(i,:)  , LS_iter(i,:)  , LS_time(i,:) ] = Matrix_Vector_Production_invQp( vec_x, 'normal', mtx.mtx_palindromic, funhand, funhand.fun_invPp, Info);
%         [ vec_y_t, relres_t(i,:), LS_iter_t(i,:), LS_time_t(i,:) ] = Matrix_Vector_Production_invQp( vec_x, 'trans' , mtx.mtx_palindromic, funhand, funhand.fun_invPp, Info);
    end
%    if strcmp(Info.linear_system.solver,'gmres')
%        LS_iter(:,3)   = Info.linear_system.restartSize*LS_iter(:,1)   + LS_iter(:,2)  ;
%         LS_iter_t(:,3) = Info.linear_system.restartSize*LS_iter_t(:,1) + LS_iter_t(:,2);
%    end
    Info.freq                 = freq;
    Info.theta_array          = theta;
    Info.ratio_array          = ratio;
    Info.linear_system.it     = LS_iter;
    Info.linear_system.relres = relres;
    Info.linear_system.time   = LS_time;
    
%     Info.linear_system.it_t     = LS_iter_t;
%     Info.linear_system.relres_t = relres_t;
%     Info.linear_system.time_t   = LS_time_t;
    
%     file_name = fullfile(FAME_option.dir_name,['result_precondition_over_theta_k',num2str(FAME_option.unit_wave_vec_idx),'_w',num2str(FAME_option.freq_idx)]);
    file_name = fullfile(FAME_option.dir_name,['result_precondition_n',num2str(N^(1/3)),'_k',num2str(FAME_option.unit_wave_vec_idx),'_w',num2str(FAME_option.freq_idx)]);
    save(file_name,'Info');
%     figure(1);
%     subplot(2,2,1)
%     surf(Theta,Ratio,reshape(relres,size(Theta)))
%     title('relres of Q_p')
%     subplot(2,2,2)
%     surf(Theta,Ratio,reshape(relres_t,size(Theta)))
%     title('relres of Q_p^T')
%     subplot(2,2,3)
%     surf(Theta,Ratio,reshape(absres,size(Theta)))
%     title('absres of Q_p')
%     subplot(2,2,4)
%     surf(Theta,Ratio,reshape(absres_t,size(Theta)))
%     title('absres of Q_p^T')
%     
%     figure(2);
%     subplot(1,2,1)
%     surf(Theta,Ratio,reshape(LS_iter(:,3),size(Theta)))
%     title('LS_iter of Q_p')
%     subplot(1,2,2)
%     surf(Theta,Ratio,reshape(LS_iter_t(:,3),size(Theta)))
%     title('LS_iter of Q_p^T')
    
    
end