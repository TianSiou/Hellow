function B = FAME_Matrix_B_Biisotropic_fem_2( B, C_1, C_2, C_3, delta_1, delta_2, delta_3 )
    N = length(C_1);
    I = speye(N);

    eps_1 = 0.5*  ( B.B_eps + B.B_eps_x);
    eps_2 = 0.5*  ( B.B_eps + B.B_eps_y);
    eps_3 = 0.5*  ( B.B_eps + B.B_eps_z);
    mu_1   = 0.25* ( B.B_mu + B.B_mu_y + B.B_mu_z + B.B_mu_yz);
    mu_2   = 0.25* ( B.B_mu + B.B_mu_z + B.B_mu_x + B.B_mu_zx);
    mu_3   = 0.25* ( B.B_mu + B.B_mu_x + B.B_mu_y + B.B_mu_xy);
    zeta_1 = 0.125*( B.B_zeta + B.B_zeta_y*(I+delta_2*C_2) + B.B_zeta_z*(I+delta_3*C_3) + B.B_zeta_yz*(I+delta_2*C_2)*(I+delta_3*C_3) ) * (2*I + delta_1*C_1');
    zeta_2 = 0.125*( B.B_zeta + B.B_zeta_x*(I+delta_1*C_1) + B.B_zeta_z*(I+delta_3*C_3) + B.B_zeta_zx*(I+delta_3*C_3)*(I+delta_1*C_1) ) * (2*I + delta_2*C_2');
    zeta_3 = 0.125*( B.B_zeta + B.B_zeta_y*(I+delta_2*C_2) + B.B_zeta_x*(I+delta_1*C_1) + B.B_zeta_xy*(I+delta_1*C_1)*(I+delta_2*C_2) ) * (2*I + delta_3*C_3');
    xi_1   = 0.125*( B.B_xi + B.B_xi_x*(I + delta_1*C_1) ) * (4*I + 2*delta_2*C_2' + 2*delta_3*C_3' + delta_2*delta_3*C_2'*C_3');
    xi_2   = 0.125*( B.B_xi + B.B_xi_y*(I + delta_2*C_2) ) * (4*I + 2*delta_1*C_1' + 2*delta_3*C_3' + delta_1*delta_3*C_3'*C_1');
    xi_3   = 0.125*( B.B_xi + B.B_xi_z*(I + delta_3*C_3) ) * (4*I + 2*delta_2*C_2' + 2*delta_1*C_1' + delta_1*delta_2*C_1'*C_2');
%     zeta_1 = 0.125*( B.B_zeta + B.B_zeta_y*(I+delta_2*C_2') + B.B_zeta_z*(I+delta_3*C_3') + B.B_zeta_yz*(I+delta_2*C_2')*(I+delta_3*C_3') ) * (2*I + 2*delta_1*C_1);
%     zeta_2 = 0.125*( B.B_zeta + B.B_zeta_x*(I+delta_1*C_1') + B.B_zeta_z*(I+delta_3*C_3') + B.B_zeta_zx*(I+delta_3*C_3')*(I+delta_1*C_1') ) * (2*I + 2*delta_2*C_2);
%     zeta_3 = 0.125*( B.B_zeta + B.B_zeta_y*(I+delta_2*C_2') + B.B_zeta_x*(I+delta_1*C_1') + B.B_zeta_xy*(I+delta_1*C_1')*(I+delta_2*C_2') ) * (2*I + 2*delta_3*C_3);
%     xi_1   = 0.125*( B.B_xi + B.B_xi_x*(I + delta_1*C_1') ) * (4*I + 2*delta_2*C_2 + 2*delta_3*C_3 + delta_2*delta_3*C_2*C_3);
%     xi_2   = 0.125*( B.B_xi + B.B_xi_y*(I + delta_2*C_2') ) * (4*I + 2*delta_1*C_1 + 2*delta_3*C_3 + delta_1*delta_3*C_3*C_1);
%     xi_3   = 0.125*( B.B_xi + B.B_xi_z*(I + delta_3*C_3') ) * (4*I + 2*delta_2*C_2 + 2*delta_1*C_1 + delta_1*delta_2*C_1*C_2);
%     zeta_1 = 0.125*( B.B_zeta*(Kx_backward + In) + B.B_zeta_y*(Kx_backward*Ky_forward + Ky_forward) + B.B_zeta_z*(Kx_backward*Kz_forward + Kz_forward) + B.B_zeta_yz*(Kx_backward*Ky_forward*Kz_forward + Ky_forward*Kz_forward) );
%     zeta_2 = 0.125*( B.B_zeta*(Ky_backward + In) + B.B_zeta_z*(Ky_backward*Kz_forward + Kz_forward) + B.B_zeta_x*(Ky_backward*Kx_forward + Kx_forward) + B.B_zeta_zx*(Ky_backward*Kz_forward*Kx_forward + Kz_forward*Kx_forward) );
%     zeta_3 = 0.125*( B.B_zeta*(Kz_backward + In) + B.B_zeta_x*(Kz_backward*Kx_forward + Kx_forward) + B.B_zeta_y*(Kz_backward*Ky_forward + Ky_forward) + B.B_zeta_xy*(Kz_backward*Kx_forward*Ky_forward + Kx_forward*Ky_forward) );
%     xi_1   = 0.125*( B.B_xi*(In + Ky_backward + Kz_backward + Ky_backward*Kz_backward) + B.B_xi_x*Kx_forward*(In + Ky_backward + Kz_backward + Ky_backward*Kz_backward) );
%     xi_2   = 0.125*( B.B_xi*(In + Kz_backward + Kx_backward + Kz_backward*Kx_backward) + B.B_xi_y*Ky_forward*(In + Kz_backward + Kx_backward + Kz_backward*Kx_backward) );
%     xi_3   = 0.125*( B.B_xi*(In + Kx_backward + Ky_backward + Kx_backward*Ky_backward) + B.B_xi_z*Kz_forward*(In + Kx_backward + Ky_backward + Kx_backward*Ky_backward) );
    
    
    B.Eps   = blkdiag(eps_1,eps_2,eps_3);
    B.Mu    = blkdiag(mu_1,mu_2,mu_3);
    B.Zeta  = blkdiag(zeta_1,zeta_2,zeta_3);
    B.Xi    = blkdiag(xi_1,xi_2,xi_3);
    B.Phi_1 = eps_1 - xi_1*zeta_1;
    B.Phi_2 = eps_2 - xi_2*zeta_2;
    B.Phi_3 = eps_3 - xi_3*zeta_3;
    B.Phi   = blkdiag(B.Phi_1,B.Phi_2,B.Phi_3);
    
%     B.B = [B.Xi, B.Mu;-B.Eps, -B.Zeta];
%     test = norm(B.Zeta - B.Xi','fro')
%     test = norm(zeta_1 - xi_1','fro')
%     test = norm(zeta_2 - xi_2','fro')
%     test = norm(zeta_3 - xi_3','fro')
%     
%     figure
%     spy(B.Zeta - B.Mu')
    
    B.lssvr = 'lu_amd';
    switch B.lssvr
        case 'lu'
    %         [ B.L, B.U ] = lu(B.Phi); 
            [ B.L_1, B.U_1 ] = lu(B.Phi_1);
            [ B.L_2, B.U_2 ] = lu(B.Phi_2);
            [ B.L_3, B.U_3 ] = lu(B.Phi_3);
        case 'chol'
    %         B.U = chol(B.Phi);
    %         B.L = B.U';
            B.U_1 = chol(B.Phi_1);
            B.U_2 = chol(B.Phi_2);
            B.U_3 = chol(B.Phi_3);
            B.L_1 = B.U_1';
            B.L_2 = B.U_2';
            B.L_3 = B.U_3';
        case 'lu_amd'
    %         B.P         = amd(B.Phi);
    %         B.invP(B.P) = 1:3*N;
    %         [ B.L, B.U ] = lu(B.Phi(B.P,B.P));  
            B.P_1         = amd(B.Phi_1);
            B.P_2         = amd(B.Phi_2);
            B.P_3         = amd(B.Phi_3);
            B.invP_1(B.P_1) = 1:N;
            B.invP_2(B.P_2) = 1:N;
            B.invP_3(B.P_3) = 1:N;
            [ B.L_1, B.U_1 ] = lu(B.Phi_1(B.P_1,B.P_1));  
            [ B.L_2, B.U_2 ] = lu(B.Phi_2(B.P_2,B.P_2));  
            [ B.L_3, B.U_3 ] = lu(B.Phi_3(B.P_3,B.P_3));  
        case 'chol_amd'
    %         B.P         = amd(B.Phi);
    %         B.invP(B.P) = 1:3*N;
    %         B.L             = chol(B.Phi(B.P,B.P),'lower');
    %         B.U             = B.L';
            B.P_1         = amd(B.Phi_1);
            B.P_2         = amd(B.Phi_2);
            B.P_3         = amd(B.Phi_3);
            B.invP_1(B.P_1) = 1:N;
            B.invP_2(B.P_2) = 1:N;
            B.invP_3(B.P_3) = 1:N;
            B.U_1         = chol(B.Phi_1(B.P_1,B.P_1));
            B.U_2         = chol(B.Phi_2(B.P_2,B.P_2));
            B.U_3         = chol(B.Phi_3(B.P_3,B.P_3));
            B.L_1         = B.U_1';
            B.L_2         = B.U_2';
            B.L_3         = B.U_3';
    end
% switch storage_format
%     case 'lu'
%         B.Phi    = B.eps - B.Xi*B.Zeta;
%         [ B.L, B.U ] = lu(B.Phi);  
%     case 'chol'
%         B.Phi    = B.eps - B.Xi*B.Zeta;
%         B.U = chol(B.Phi);
%         B.L = B.U';
%     case 'lu_amd'
%         B.Phi       = B.eps - B.Xi*B.Zeta;
%         B.P         = amd(B.Phi);
%         B.invP(B.P) = 1:3*N;
%         [ B.L, B.U ] = lu(B.Phi(B.P,B.P));  
%     case 'chol_amd'
%         B.Phi       = B.eps - B.Xi*B.Zeta;
%         B.P         = amd(B.Phi);
%         B.invP(B.P) = 1:3*N;
%         B.L             = chol(B.Phi(B.P,B.P),'lower');
%         B.U             = B.L';
%     otherwise
%         B.Phi    = B.eps - B.Xi*B.Zeta;
% end

% Plot sparsity of B's
% h = figure(3); cla
% subplot(2,2,1); hold on
% plot([  N+1/2,    N+1/2],[      0,  3*N+1/2],'k-')
% plot([2*N+1/2,  2*N+1/2],[      0,  3*N+1/2],'k-')
% plot([      0,  3*N+1/2],[  N+1/2,    N+1/2],'k-')
% plot([      0,  3*N+1/2],[2*N+1/2,  2*N+1/2],'k-')
% plot([      0,  3*N+1/2],[      0,        0],'k-')
% plot([3*N+1/2,  3*N+1/2],[      0,  3*N+1/2],'k-')
% spy(B.Zeta)
% title('Nonzeros of \zeta_d')
% subplot(2,2,2); hold on
% plot([  N+1/2,    N+1/2],[      0,  3*N+1/2],'k-')
% plot([2*N+1/2,  2*N+1/2],[      0,  3*N+1/2],'k-')
% plot([      0,  3*N+1/2],[  N+1/2,    N+1/2],'k-')
% plot([      0,  3*N+1/2],[2*N+1/2,  2*N+1/2],'k-')
% plot([      0,  3*N+1/2],[      0,        0],'k-')
% plot([3*N+1/2,  3*N+1/2],[      0,  3*N+1/2],'k-')
% spy(B.Mu)
% title('Nonzeros of \mu_d')
% subplot(2,2,3); hold on
% plot([  N+1/2,    N+1/2],[      0,  3*N+1/2],'k-')
% plot([2*N+1/2,  2*N+1/2],[      0,  3*N+1/2],'k-')
% plot([      0,  3*N+1/2],[  N+1/2,    N+1/2],'k-')
% plot([      0,  3*N+1/2],[2*N+1/2,  2*N+1/2],'k-')
% plot([      0,  3*N+1/2],[      0,        0],'k-')
% plot([3*N+1/2,  3*N+1/2],[      0,  3*N+1/2],'k-')
% spy(B.eps)
% title('Nonzeros of \epsilon_d')
% subplot(2,2,4); hold on
% plot([  N+1/2,    N+1/2],[      0,  3*N+1/2],'k-')
% plot([2*N+1/2,  2*N+1/2],[      0,  3*N+1/2],'k-')
% plot([      0,  3*N+1/2],[  N+1/2,    N+1/2],'k-')
% plot([      0,  3*N+1/2],[2*N+1/2,  2*N+1/2],'k-')
% plot([      0,  3*N+1/2],[      0,        0],'k-')
% plot([3*N+1/2,  3*N+1/2],[      0,  3*N+1/2],'k-')
% spy(B.Xi)
% title('Nonzeros of \xi_d')
% saveas(h,['Sparsity_of_B_[',num2str(wave_vec(1)),',',num2str(wave_vec(2)),',',num2str(wave_vec(3)),'].fig'])
% close(h)
% %% Plot sparsity of Phi
% h = figure(4); cla
% % subplot(1,2,1); 
% hold on
% plot([      0,  3*N+1/2],[      0,        0],'k-')
% plot([  N+1/2,    N+1/2],[      0,  3*N+1/2],'k-')
% plot([2*N+1/2,  2*N+1/2],[      0,  3*N+1/2],'k-')
% plot([      0,  3*N+1/2],[  N+1/2,    N+1/2],'k-')
% plot([      0,  3*N+1/2],[2*N+1/2,  2*N+1/2],'k-')
% plot([3*N+1/2,  3*N+1/2],[      0,  3*N+1/2],'k-')
% spy(B.Phi)
% title('Nonzeros of \Phi_d')
% figure(5); cla
% % subplot(1,2,2); 
% B.U_11 = chol(B.Phi_1);
% B.U_22 = chol(B.Phi_2);
% B.U_33 = chol(B.Phi_3);
% hold on
% plot([      0,  3*N+1/2],[      0,        0],'k-')
% plot([  N+1/2,    N+1/2],[      0,  3*N+1/2],'k-')
% plot([2*N+1/2,  2*N+1/2],[      0,  3*N+1/2],'k-')
% plot([      0,  3*N+1/2],[  N+1/2,    N+1/2],'k-')
% plot([      0,  3*N+1/2],[2*N+1/2,  2*N+1/2],'k-')
% plot([3*N+1/2,  3*N+1/2],[      0,  3*N+1/2],'k-')
% spy(blkdiag(B.U_11,B.U_22,B.U_33))
% title('Nonzeros of U without amd')
% figure(6); cla
% % subplot(1,2,2); 
% hold on
% plot([      0,  3*N+1/2],[      0,        0],'k-')
% plot([  N+1/2,    N+1/2],[      0,  3*N+1/2],'k-')
% plot([2*N+1/2,  2*N+1/2],[      0,  3*N+1/2],'k-')
% plot([      0,  3*N+1/2],[  N+1/2,    N+1/2],'k-')
% plot([      0,  3*N+1/2],[2*N+1/2,  2*N+1/2],'k-')
% plot([3*N+1/2,  3*N+1/2],[      0,  3*N+1/2],'k-')
% spy(blkdiag(B.U_1,B.U_2,B.U_3))
% title('Nonzeros of U with amd')
% saveas(h,['Sparsity_of_Phi_[',num2str(wave_vec(1)),',',num2str(wave_vec(2)),',',num2str(wave_vec(3)),'].fig'])
% close(h)
%% Check hermitian of xi and zeta
% test_hermitian_xi_zeta = norm(B.Zeta - B.Xi','fro')
%% Test positive definite of Phi
% Phi = Eps - Zeta*Xi;
% eig(full(Phi))

end