function Plot_Small_Dimension(C, C_1, C_2, C_3, K_tilde, K_1_tilde, K_2_tilde, K_3_tilde, B, Lambdas, freq, N, Par_mesh, Par_lattice)
    
    % 選擇 Symplectic GEP 的 shift 值 \eta_0
    Info.shift_ss = -1.8;
    % 計算 PQEP 的 shift 值 nu_0
    Info.shift_palindromic =  0.5*(Info.shift_ss + sqrt(Info.shift_ss^2 - 4)); % nu_0
    % 調整內部線性系統的選項
    Info.linear_system.dimension   = length(C);
    Info.linear_system.tol         = 1e-14;
    Info.linear_system.solver      = 'gmres';
    Info.linear_system.v0          = rand(3*N,1);
    Info.linear_system.isprecond   = 1;
    Info.linear_system.maxit       = 100;
    Info.linear_system.restartSize = 30;
    % 選擇線性系統 Qp*x = b 的 precondition 參數 a_eps
    Info.linear_system.alpha = 0.25*(Info.shift_palindromic-1)^2*freq^2*sum(B.B_eps)/length(B.B_eps);
%     option.alpha = 1000*option.alpha;
    % 比較不同 a_eps 的效果
%     NumericalResult_compare_precondition
    
    % Construst the coefficient matrices of GQEP, PQEP, and their transformations
    mtx = Matrix_Generate_1(C, C_1, C_2, C_3, K_tilde, K_1_tilde, K_2_tilde, K_3_tilde, B, Lambdas, N);
    mtx = Matrix_Generate_2(mtx, freq, Info.shift_palindromic, Info.linear_system.alpha );    
    % 建構各個 function handle
    funhand = Funhand_Generate(C, K_tilde,B, mtx, Par_mesh, Par_lattice, Lambdas, freq, Info, N);
    % 在小矩陣時畫出各種轉化下的 Spectrum
    Test_Dense_Eigenproblems(mtx.mtx_gyroscopic, mtx.mtx_SHH, mtx.mtx_palindromic, mtx.mtx_symplecitc, mtx.mtx_ss, mtx.mtx_SH);       
end