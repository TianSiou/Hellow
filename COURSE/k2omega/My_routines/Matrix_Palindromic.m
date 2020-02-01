function mtx_palindromic = Matrix_Palindromic( C, C_1, C_2, C_3, K_tilde, K_1_tilde, K_2_tilde, K_3_tilde, mtx_gyroscopic, Lambdas )
    N = length(Lambdas.Lambda_x);

    mtx_palindromic.fun_A = @(freq) 0.25*(mtx_gyroscopic.fun_M(freq)-mtx_gyroscopic.G+mtx_gyroscopic.K);
    mtx_palindromic.fun_Q = @(freq) 0.5 *(mtx_gyroscopic.fun_M(freq)-mtx_gyroscopic.K);
    
%     mtx_palindromic.fun_Qp = @(sigma,freq) sigma^2*mtx_palindromic.fun_A(freq).'- sigma*mtx_palindromic.fun_Q(freq) + mtx_palindromic.fun_A(freq);
%     mtx_palindromic.fun_Pp = @(sigma,alpha) 0.25*(sigma^2-2*sigma+1)*(C.'*C) + 0.25*(sigma^2-1)*mtx_gyroscopic.G + 0.25*(sigma^2+2*sigma+1)*mtx_gyroscopic.K - alpha*speye(3*N);
%     mtx_palindromic.Pp_funhand  = @(vec_x, sigma,alpha) 0.25*(sigma^2-2*sigma+1)*(C.'*(C*vec_x)) + 0.25*(sigma^2-1)*(mtx_gyroscopic.G*vec_x) + 0.25*(sigma^2+2*sigma+1)*mtx_gyroscopic.K*vec_x - alpha*vec_x;
%     mtx_palindromic.Ppt_funhand = @(vec_x, sigma,alpha) 0.25*(sigma^2-2*sigma+1)*(C.'*(C*vec_x)) - 0.25*(sigma^2-1)*(mtx_gyroscopic.G*vec_x) + 0.25*(sigma^2+2*sigma+1)*mtx_gyroscopic.K*vec_x - alpha*vec_x;
    
    Lambda_x = Lambdas.Lambda_x;
    Lambda_y = Lambdas.Lambda_y;
    Lambda_z = Lambdas.Lambda_z;
    Lambda_x_tilde = Lambdas.Lambda_x_tilde;
    Lambda_y_tilde = Lambdas.Lambda_y_tilde;
    Lambda_z_tilde = Lambdas.Lambda_z_tilde;   
    
%     O_N = sparse(N,N);
%     Lambda   = [             O_N, -spdiags(Lambda_z,0,N,N),  spdiags(Lambda_y,0,N,N);
%                   spdiags(Lambda_z,0,N,N),             O_N, -spdiags(Lambda_x,0,N,N);
%                  -spdiags(Lambda_y,0,N,N),  spdiags(Lambda_x,0,N,N),             O_N];
%     Lambda_tilde   = [             O_N, -spdiags(Lambda_z_tilde,0,N,N),  spdiags(Lambda_y_tilde,0,N,N);
%                         spdiags(Lambda_z_tilde,0,N,N),             O_N, -spdiags(Lambda_x_tilde,0,N,N);
%                        -spdiags(Lambda_y_tilde,0,N,N),  spdiags(Lambda_x_tilde,0,N,N),             O_N];
%     N1 = @(sigma) 0.5*(sigma-1)*[C_1;C_2;C_3] + 0.5*(sigma+1)*[K_1_tilde;K_2_tilde;K_3_tilde];  
%     N2 = @(sigma) 0.5*(conj(sigma)-1)*[C_1;C_2;C_3] - 0.5*(conj(sigma)+1)*[K_1_tilde;K_2_tilde;K_3_tilde];
    Lambda_N1_x = @(sigma) 0.5*(sigma-1)*Lambda_x + 0.5*(sigma+1)*Lambda_x_tilde;  Lambda_N2_x = @(sigma) 0.5*(conj(sigma)-1)*Lambda_x - 0.5*(conj(sigma)+1)*Lambda_x_tilde;
    Lambda_N1_y = @(sigma) 0.5*(sigma-1)*Lambda_y + 0.5*(sigma+1)*Lambda_y_tilde;  Lambda_N2_y = @(sigma) 0.5*(conj(sigma)-1)*Lambda_y - 0.5*(conj(sigma)+1)*Lambda_y_tilde;
    Lambda_N1_z = @(sigma) 0.5*(sigma-1)*Lambda_z + 0.5*(sigma+1)*Lambda_z_tilde;  Lambda_N2_z = @(sigma) 0.5*(conj(sigma)-1)*Lambda_z - 0.5*(conj(sigma)+1)*Lambda_z_tilde;
    Lambda_N1bar_x = @(sigma) 0.5*(conj(sigma)-1)*Lambda_x + 0.5*(conj(sigma)+1)*Lambda_x_tilde;  Lambda_N2bar_x = @(sigma) 0.5*(sigma-1)*Lambda_x - 0.5*(sigma+1)*Lambda_x_tilde;
    Lambda_N1bar_y = @(sigma) 0.5*(conj(sigma)-1)*Lambda_y + 0.5*(conj(sigma)+1)*Lambda_y_tilde;  Lambda_N2bar_y = @(sigma) 0.5*(sigma-1)*Lambda_y - 0.5*(sigma+1)*Lambda_y_tilde;
    Lambda_N1bar_z = @(sigma) 0.5*(conj(sigma)-1)*Lambda_z + 0.5*(conj(sigma)+1)*Lambda_z_tilde;  Lambda_N2bar_z = @(sigma) 0.5*(sigma-1)*Lambda_z - 0.5*(sigma+1)*Lambda_z_tilde;
    % test nu_0 = -1
%     mtx_palindromic.Nc               = [C_1;C_2;C_3];
%     mtx_palindromic.Lambda_q         = [0;Lambdas.Lambda_q];
%     mtx_palindromic.Lambda_c         = [spdiags(Lambdas.Lambda_x,0,N,N);spdiags(Lambdas.Lambda_y,0,N,N);spdiags(Lambdas.Lambda_z,0,N,N)];
%     mtx_palindromic.fun_Sigma_pp     = @(alpha) ( [mtx_palindromic.Lambda_q;mtx_palindromic.Lambda_q;mtx_palindromic.Lambda_q] - alpha );
%     mtx_palindromic.fun_invSigma_pp  = @(alpha) 1./( [mtx_palindromic.Lambda_q;mtx_palindromic.Lambda_q;mtx_palindromic.Lambda_q] - alpha );
    
    Lambda_q     = @(sigma) conj(Lambda_N2_x(sigma)).*Lambda_N1_x(sigma) + conj(Lambda_N2_y(sigma)).*Lambda_N1_y(sigma) + conj(Lambda_N2_z(sigma)).*Lambda_N1_z(sigma);
    Lambda_qt    = @(sigma) conj(Lambda_N1bar_x(sigma)).*Lambda_N2bar_x(sigma) + conj(Lambda_N1bar_y(sigma)).*Lambda_N2bar_y(sigma) + conj(Lambda_N1bar_z(sigma)).*Lambda_N2bar_z(sigma);
    mtx_palindromic.fun_Lambda_N1  = @(sigma) [spdiags(Lambda_N1_x(sigma),0,N,N);spdiags(Lambda_N1_y(sigma),0,N,N);spdiags(Lambda_N1_z(sigma),0,N,N)];
    mtx_palindromic.fun_Lambda_N2  = @(sigma) [spdiags(Lambda_N2_x(sigma),0,N,N);spdiags(Lambda_N2_y(sigma),0,N,N);spdiags(Lambda_N2_z(sigma),0,N,N)];
    mtx_palindromic.fun_Lambda_N1bar  = @(sigma) [spdiags(Lambda_N1bar_x(sigma),0,N,N);spdiags(Lambda_N1bar_y(sigma),0,N,N);spdiags(Lambda_N1bar_z(sigma),0,N,N)];
    mtx_palindromic.fun_Lambda_N2bar  = @(sigma) [spdiags(Lambda_N2bar_x(sigma),0,N,N);spdiags(Lambda_N2bar_y(sigma),0,N,N);spdiags(Lambda_N2bar_z(sigma),0,N,N)];
%     mtx_palindromic.fun_Sigma_p  = @(sigma,alpha) ( [Lambda_q(sigma);Lambda_q(sigma);Lambda_q(sigma)] - alpha );
%     mtx_palindromic.fun_Sigma_pt = @(sigma,alpha) ( [Lambda_qt(sigma);Lambda_qt(sigma);Lambda_qt(sigma)] - alpha );
    mtx_palindromic.fun_invSigma_p  = @(sigma,alpha) 1./( [Lambda_q(sigma);Lambda_q(sigma);Lambda_q(sigma)] - alpha );
    mtx_palindromic.fun_invSigma_pt = @(sigma,alpha) 1./( [Lambda_qt(sigma);Lambda_qt(sigma);Lambda_qt(sigma)] - alpha );
%     mtx_palindromic.fun_D  = @(sigma) ( 0.5*(conj(sigma)-1)*Lambda - 0.5*(conj(sigma)+1)*Lambda_tilde )'*( 0.5*(sigma-1)*Lambda + 0.5*(sigma+1)*Lambda_tilde);
%     test = -2;
%     alpha = 3;
%     norm(mtx_palindromic.fun_Pp(test,alpha) - blkdiag(N2(test)'*N1(test),N2(test)'*N1(test),N2(test)'*N1(test)) ...
%         + N1(test)*N2(test)'+alpha*speye(3*N),'fro')
%     norm(mtx_palindromic.fun_D(test) - blkdiag(mtx_palindromic.fun_Lambda_N2(test)'*mtx_palindromic.fun_Lambda_N1(test),mtx_palindromic.fun_Lambda_N2(test)'*mtx_palindromic.fun_Lambda_N1(test),mtx_palindromic.fun_Lambda_N2(test)'*mtx_palindromic.fun_Lambda_N1(test)) ...
%         + mtx_palindromic.fun_Lambda_N1(test)*mtx_palindromic.fun_Lambda_N2(test)','fro')
    
%     Dp_0.d11 = -dot(Lambda_z_tilde,Lambda_z_tilde,2) - dot(Lambda_y_tilde,Lambda_y_tilde,2);
%     Dp_0.d12 =  dot(Lambda_y_tilde,Lambda_x_tilde,2);
%     Dp_0.d13 =  dot(Lambda_z_tilde,Lambda_x_tilde,2);
%     Dp_0.d22 = -dot(Lambda_z_tilde,Lambda_z_tilde,2) - dot(Lambda_x_tilde,Lambda_x_tilde,2);
%     Dp_0.d23 =  dot(Lambda_z_tilde,Lambda_y_tilde,2);
%     Dp_0.d33 = -dot(Lambda_y_tilde,Lambda_y_tilde,2) - dot(Lambda_x_tilde,Lambda_x_tilde,2);
% 
%     Dp_1.d11 = -(dot(Lambda_z_tilde,Lambda_z,2) + dot(Lambda_y_tilde,Lambda_y,2)) + ...
%                 (dot(Lambda_z,Lambda_z_tilde,2) + dot(Lambda_y,Lambda_y_tilde,2)) ;
%     Dp_1.d12 =  dot(Lambda_y_tilde,Lambda_x,2) - dot(Lambda_y,Lambda_x_tilde,2);
%     Dp_1.d13 =  dot(Lambda_z_tilde,Lambda_x,2) - dot(Lambda_z,Lambda_x_tilde,2);
%     Dp_1.d22 = -(dot(Lambda_z_tilde,Lambda_z,2) + dot(Lambda_x_tilde,Lambda_x,2)) + ...
%                 (dot(Lambda_z,Lambda_z_tilde,2) + dot(Lambda_x,Lambda_x_tilde,2)) ;
%     Dp_1.d23 =  dot(Lambda_z_tilde,Lambda_y,2) - dot(Lambda_z,Lambda_y_tilde,2);
%     Dp_1.d33 = -(dot(Lambda_y_tilde,Lambda_y,2) + dot(Lambda_x_tilde,Lambda_x,2)) + ...
%                 (dot(Lambda_y,Lambda_y_tilde,2) + dot(Lambda_x,Lambda_x_tilde,2)) ;
%     
%     Dp_2.d11_fun = @(freq,alpha)  dot(Lambda_z,Lambda_z,2) + dot(Lambda_y,Lambda_y,2) - alpha;
%     Dp_2.d12     = -dot(Lambda_y,Lambda_x,2);
%     Dp_2.d13     = -dot(Lambda_z,Lambda_x,2);
%     Dp_2.d22_fun = @(freq,alpha)  dot(Lambda_z,Lambda_z,2) + dot(Lambda_x,Lambda_x,2) - alpha;
%     Dp_2.d23     = -dot(Lambda_z,Lambda_y,2);
%     Dp_2.d33_fun = @(freq,alpha)  dot(Lambda_y,Lambda_y,2) + dot(Lambda_x,Lambda_x,2) - alpha;
    % function handles for Dp
%     mtx_palindromic.fun_Dp.fun_d11 = @(sigma,freq,alpha) 0.25*(sigma^2-2*sigma+1)*Dp_2.d11_fun(freq,alpha) + ...
%                                                          0.25*(sigma^2-1)        *Dp_1.d11 + ...
%                                                          0.25*(sigma^2+2*sigma+1)*Dp_0.d11;
%     mtx_palindromic.fun_Dp.fun_d12 = @(sigma) 0.25*(sigma^2-2*sigma+1)*Dp_2.d12 + ...
%                                               0.25*(sigma^2-1)        *Dp_1.d12 + ...
%                                               0.25*(sigma^2+2*sigma+1)*Dp_0.d12;
%     mtx_palindromic.fun_Dp.fun_d13 = @(sigma) 0.25*(sigma^2-2*sigma+1)*Dp_2.d13 + ...
%                                               0.25*(sigma^2-1)        *Dp_1.d13 + ...
%                                               0.25*(sigma^2+2*sigma+1)*Dp_0.d13;
%     mtx_palindromic.fun_Dp.fun_d21 = @(sigma) 0.25*(sigma^2-2*sigma+1)*conj(Dp_2.d12) + ...
%                                               0.25*(sigma^2-1)        *(-conj(Dp_1.d12)) + ...
%                                               0.25*(sigma^2+2*sigma+1)*conj(Dp_0.d12);
%     mtx_palindromic.fun_Dp.fun_d22 = @(sigma,freq,alpha) 0.25*(sigma^2-2*sigma+1)*Dp_2.d22_fun(freq,alpha) + ...
%                                                          0.25*(sigma^2-1)*Dp_1.d22 + ...
%                                                          0.25*(sigma^2+2*sigma+1)*Dp_0.d22;
%     mtx_palindromic.fun_Dp.fun_d23 = @(sigma) 0.25*(sigma^2-2*sigma+1)*Dp_2.d23 + ...
%                                               0.25*(sigma^2-1)        *Dp_1.d23 + ...
%                                               0.25*(sigma^2+2*sigma+1)*Dp_0.d23;
%     mtx_palindromic.fun_Dp.fun_d31 = @(sigma) 0.25*(sigma^2-2*sigma+1)*conj(Dp_2.d13) + ...
%                                               0.25*(sigma^2-1)        *(-conj(Dp_1.d13)) + ...
%                                               0.25*(sigma^2+2*sigma+1)*conj(Dp_0.d13);
%     mtx_palindromic.fun_Dp.fun_d32 = @(sigma) 0.25*(sigma^2-2*sigma+1)*conj(Dp_2.d23) + ...
%                                               0.25*(sigma^2-1)        *(-conj(Dp_1.d23)) + ...
%                                               0.25*(sigma^2+2*sigma+1)*conj(Dp_0.d23);
%     mtx_palindromic.fun_Dp.fun_d33 = @(sigma,freq,alpha) 0.25*(sigma^2-2*sigma+1)*Dp_2.d33_fun(freq,alpha) + ...
%                                                          0.25*(sigma^2-1)*Dp_1.d33 + ...
%                                                          0.25*(sigma^2+2*sigma+1)*Dp_0.d33;
%     mtx_palindromic.fun_Dp.fun_Dp = @(sigma,freq,alpha) ...
%         [spdiags(mtx_palindromic.fun_Dp.fun_d11(sigma,freq,alpha),0,N,N), spdiags(mtx_palindromic.fun_Dp.fun_d12(sigma),0,N,N), spdiags(mtx_palindromic.fun_Dp.fun_d13(sigma),0,N,N);
%          spdiags(mtx_palindromic.fun_Dp.fun_d21(sigma),0,N,N), spdiags(mtx_palindromic.fun_Dp.fun_d22(sigma,freq,alpha),0,N,N), spdiags(mtx_palindromic.fun_Dp.fun_d23(sigma),0,N,N);
%          spdiags(mtx_palindromic.fun_Dp.fun_d31(sigma),0,N,N), spdiags(mtx_palindromic.fun_Dp.fun_d32(sigma),0,N,N), spdiags(mtx_palindromic.fun_Dp.fun_d33(sigma,freq,alpha),0,N,N)];
    % function handles for DpT
%     mtx_palindromic.fun_DpT.fun_d11 = @(sigma,freq,alpha) 0.25*(sigma^2-2*sigma+1)*Dp_2.d11_fun(freq,alpha) + ...
%                                                          0.25*(sigma^2-1)        *(-Dp_1.d11) + ...
%                                                          0.25*(sigma^2+2*sigma+1)*Dp_0.d11;
%     mtx_palindromic.fun_DpT.fun_d12 = @(sigma) 0.25*(sigma^2-2*sigma+1)*Dp_2.d12 + ...
%                                               0.25*(sigma^2-1)        *(-Dp_1.d12) + ...
%                                               0.25*(sigma^2+2*sigma+1)*Dp_0.d12;
%     mtx_palindromic.fun_DpT.fun_d13 = @(sigma) 0.25*(sigma^2-2*sigma+1)*Dp_2.d13 + ...
%                                               0.25*(sigma^2-1)        *(-Dp_1.d13) + ...
%                                               0.25*(sigma^2+2*sigma+1)*Dp_0.d13;
%     mtx_palindromic.fun_DpT.fun_d21 = @(sigma) 0.25*(sigma^2-2*sigma+1)*conj(Dp_2.d12) + ...
%                                               0.25*(sigma^2-1)        *(conj(Dp_1.d12)) + ...
%                                               0.25*(sigma^2+2*sigma+1)*conj(Dp_0.d12);
%     mtx_palindromic.fun_DpT.fun_d22 = @(sigma,freq,alpha) 0.25*(sigma^2-2*sigma+1)*Dp_2.d22_fun(freq,alpha) + ...
%                                                          0.25*(sigma^2-1)*(-Dp_1.d22) + ...
%                                                          0.25*(sigma^2+2*sigma+1)*Dp_0.d22;
%     mtx_palindromic.fun_DpT.fun_d23 = @(sigma) 0.25*(sigma^2-2*sigma+1)*Dp_2.d23 + ...
%                                               0.25*(sigma^2-1)        *(-Dp_1.d23) + ...
%                                               0.25*(sigma^2+2*sigma+1)*Dp_0.d23;
%     mtx_palindromic.fun_DpT.fun_d31 = @(sigma) 0.25*(sigma^2-2*sigma+1)*conj(Dp_2.d13) + ...
%                                               0.25*(sigma^2-1)        *conj(Dp_1.d13) + ...
%                                               0.25*(sigma^2+2*sigma+1)*conj(Dp_0.d13);
%     mtx_palindromic.fun_DpT.fun_d32 = @(sigma) 0.25*(sigma^2-2*sigma+1)*conj(Dp_2.d23) + ...
%                                               0.25*(sigma^2-1)        *conj(Dp_1.d23) + ...
%                                               0.25*(sigma^2+2*sigma+1)*conj(Dp_0.d23);
%     mtx_palindromic.fun_DpT.fun_d33 = @(sigma,freq,alpha) 0.25*(sigma^2-2*sigma+1)*Dp_2.d33_fun(freq,alpha) + ...
%                                                          0.25*(sigma^2-1)*(-Dp_1.d33) + ...
%                                                          0.25*(sigma^2+2*sigma+1)*Dp_0.d33;
%     mtx_palindromic.fun_DpT.fun_DpT = @(sigma,freq,alpha) ...
%         [spdiags(mtx_palindromic.fun_DpT.fun_d11(sigma,freq,alpha),0,N,N), spdiags(mtx_palindromic.fun_DpT.fun_d12(sigma),0,N,N), spdiags(mtx_palindromic.fun_DpT.fun_d13(sigma),0,N,N);
%          spdiags(mtx_palindromic.fun_DpT.fun_d21(sigma),0,N,N), spdiags(mtx_palindromic.fun_DpT.fun_d22(sigma,freq,alpha),0,N,N), spdiags(mtx_palindromic.fun_DpT.fun_d23(sigma),0,N,N);
%          spdiags(mtx_palindromic.fun_DpT.fun_d31(sigma),0,N,N), spdiags(mtx_palindromic.fun_DpT.fun_d32(sigma),0,N,N), spdiags(mtx_palindromic.fun_DpT.fun_d33(sigma,freq,alpha),0,N,N)];
%     
%     mtx_palindromic.Pi = kron(1:N,[1,1,1]) + kron(ones(1,N),[0,N,2*N]);
end

