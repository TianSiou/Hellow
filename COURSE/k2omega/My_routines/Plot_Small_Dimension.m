function Plot_Small_Dimension(C, C_1, C_2, C_3, K_tilde, K_1_tilde, K_2_tilde, K_3_tilde, B, Lambdas, freq, N, Par_mesh, Par_lattice)
    
    % ��� Symplectic GEP �� shift �� \eta_0
    Info.shift_ss = -1.8;
    % �p�� PQEP �� shift �� nu_0
    Info.shift_palindromic =  0.5*(Info.shift_ss + sqrt(Info.shift_ss^2 - 4)); % nu_0
    % �վ㤺���u�ʨt�Ϊ��ﶵ
    Info.linear_system.dimension   = length(C);
    Info.linear_system.tol         = 1e-14;
    Info.linear_system.solver      = 'gmres';
    Info.linear_system.v0          = rand(3*N,1);
    Info.linear_system.isprecond   = 1;
    Info.linear_system.maxit       = 100;
    Info.linear_system.restartSize = 30;
    % ��ܽu�ʨt�� Qp*x = b �� precondition �Ѽ� a_eps
    Info.linear_system.alpha = 0.25*(Info.shift_palindromic-1)^2*freq^2*sum(B.B_eps)/length(B.B_eps);
%     option.alpha = 1000*option.alpha;
    % ������P a_eps ���ĪG
%     NumericalResult_compare_precondition
    
    % Construst the coefficient matrices of GQEP, PQEP, and their transformations
    mtx = Matrix_Generate_1(C, C_1, C_2, C_3, K_tilde, K_1_tilde, K_2_tilde, K_3_tilde, B, Lambdas, N);
    mtx = Matrix_Generate_2(mtx, freq, Info.shift_palindromic, Info.linear_system.alpha );    
    % �غc�U�� function handle
    funhand = Funhand_Generate(C, K_tilde,B, mtx, Par_mesh, Par_lattice, Lambdas, freq, Info, N);
    % �b�p�x�}�ɵe�X�U����ƤU�� Spectrum
    Test_Dense_Eigenproblems(mtx.mtx_gyroscopic, mtx.mtx_SHH, mtx.mtx_palindromic, mtx.mtx_symplecitc, mtx.mtx_ss, mtx.mtx_SH);       
end