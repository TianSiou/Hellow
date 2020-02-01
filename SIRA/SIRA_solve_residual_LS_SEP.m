function [ vec_t ] = SIRA_solve_residual_LS_SEP(A, sigma, rhs, LS_solver, pcg_iter_number, ...
    linear_system_tol, Sigma_r,init_vec, prec_M)
%Solve the linear system
%        ( A - tau I ) vec_t = rhs.
%    

% global inner_iter_count outer_iter_count
% fun_shift_Ar   = @(vec_x) A(vec_x) - sigma * vec_x;
% global inner_iter inner_count inner_cpu_time 
% inner_count = inner_count + 1;
if ( nargin <= 6 )
    if isa(A,'function_handle')
        switch LS_solver
            case 'pcg'
%                 [vec_t, flag, RELRES, iter_num ] = pcg(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
%                      linear_system_tol, pcg_iter_number); 
                [vec_t, flag, RELRES, iter_num ] = pcg(fun_shift_Ar, rhs, ...
                     linear_system_tol, pcg_iter_number); 
%                 [vec_t, flag, RELRES, iter_num ] = pcg(A,rhs, ...
%                      linear_system_tol, pcg_iter_number); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
    
            case 'minres'
%                 [vec_t, flag, RELRES, iter_num ] = minres(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
%                      linear_system_tol, pcg_iter_number); 
                [vec_t, flag, RELRES, iter_num ] = minres(fun_shift_Ar, rhs, ...
                     linear_system_tol, pcg_iter_number); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'bicgstab'
                [vec_t, flag, RELRES, iter_num ] = bicgstab(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
                     linear_system_tol, pcg_iter_number); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'bicgstabl'
                [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
                     linear_system_tol, pcg_iter_number); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'cgs'
                [vec_t, flag, RELRES, iter_num ] = cgs(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
                     linear_system_tol, pcg_iter_number); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'gmres'
                RESTRT       = 10;
                [vec_t, flag, RELRES, iter_num ] = gmres(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
                     RESTRT, linear_system_tol, pcg_iter_number); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num(1)*RESTRT+iter_num(2), RELRES, flag);
            
            otherwise
                [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
                     linear_system_tol, pcg_iter_number); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        end
    else
        mtx_A = A - sigma * speye(size(A,1)); 
        switch LS_solver
            case 'pcg' 
                [vec_t, flag, RELRES, iter_num ] = pcg(mtx_A, rhs, linear_system_tol, pcg_iter_number); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
    
            case 'minres'
                [vec_t, flag, RELRES, iter_num ] = minres(mtx_A, rhs, linear_system_tol, pcg_iter_number); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'bicgstab'
                [vec_t, flag, RELRES, iter_num ] = bicgstab(mtx_A, rhs, linear_system_tol, pcg_iter_number); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'bicgstabl'
                [vec_t, flag, RELRES, iter_num ] = bicgstabl(mtx_A, rhs, linear_system_tol, pcg_iter_number); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'cgs'
                [vec_t, flag, RELRES, iter_num ] = cgs(mtx_A, rhs, linear_system_tol, pcg_iter_number); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'gmres'
                RESTRT       = 10;
                [vec_t, flag, RELRES, iter_num ] = gmres(mtx_A, rhs, RESTRT, linear_system_tol, pcg_iter_number); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num(1)*RESTRT+iter_num(2), RELRES, flag);
            
            otherwise
                [vec_t, flag, RELRES, iter_num ] = bicgstabl(mtx_A, rhs, linear_system_tol, pcg_iter_number); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        end
    end
elseif ( nargin == 7 )
    if isa(A,'function_handle')
        switch LS_solver
            case 'pcg'
%                 [vec_t, flag, RELRES, iter_num ] = pcg(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
%                      linear_system_tol, pcg_iter_number, [], [], init_vec); 
     
                rhs  = Sigma_r .\ rhs;
%                 time_start_pcg = tic;
%                 [vec_t, flag, RELRES, inner_iter(inner_count) ] = pcg(A, rhs, ...
%                      linear_system_tol, pcg_iter_number);
%                  inner_cpu_time(inner_count) = toc(time_start_pcg);
%                  vec_t = Sigma_r .\ vec_t;
%                 fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',inner_iter(inner_count), RELRES, flag);
                 [vec_t, flag, RELRES, iter_num ] = pcg(A, rhs, ...
                                     linear_system_tol, pcg_iter_number);
%                                  inner_cpu_time(inner_count) = toc(time_start_pcg);
                                 vec_t = Sigma_r .\ vec_t;
                                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
    
            case 'minres'
%                 [vec_t, flag, RELRES, iter_num ] = minres(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
%                      linear_system_tol, pcg_iter_number, [], [], init_vec); 
                 rhs  = Sigma_r .\ rhs;
                [vec_t, flag, RELRES, iter_num ] = minres(@(vec_x) A(vec_x) - sigma* ( (Sigma_r.^2) .\ vec_x),rhs, ...
                     linear_system_tol, pcg_iter_number);
                 vec_t = Sigma_r .\ vec_t;
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'bicgstab'
                [vec_t, flag, RELRES, iter_num ] = bicgstab(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
                     linear_system_tol, pcg_iter_number, [], [], init_vec); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'bicgstabl'
                [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
                     linear_system_tol, pcg_iter_number, [], [], init_vec); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'cgs'
                [vec_t, flag, RELRES, iter_num ] = cgs(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
                     linear_system_tol, pcg_iter_number, [], [], init_vec); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'gmres'
                RESTRT       = 10;
                [vec_t, flag, RELRES, iter_num ] = gmres(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
                     RESTRT, linear_system_tol, pcg_iter_number, [], [], init_vec); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num(1)*RESTRT+iter_num(2), RELRES, flag);
            
            otherwise
                [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
                     linear_system_tol, pcg_iter_number, [], [], init_vec); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        end
    else
        mtx_A = A - sigma * speye(size(A,1));
        switch LS_solver
            case 'pcg'
                [vec_t, flag, RELRES, iter_num ] = pcg(mtx_A, rhs, linear_system_tol, pcg_iter_number, [], [], init_vec); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
    
            case 'minres'
                [vec_t, flag, RELRES, iter_num ] = minres(mtx_A, rhs, linear_system_tol, pcg_iter_number, [], [], init_vec); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'bicgstab'
                [vec_t, flag, RELRES, iter_num ] = bicgstab(mtx_A, rhs, linear_system_tol, pcg_iter_number, [], [], init_vec); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'bicgstabl'
                [vec_t, flag, RELRES, iter_num ] = bicgstabl(mtx_A, rhs, linear_system_tol, pcg_iter_number, [], [], init_vec); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'cgs'
                [vec_t, flag, RELRES, iter_num ] = cgs(mtx_A, rhs, linear_system_tol, pcg_iter_number, [], [], init_vec); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'gmres'
                RESTRT       = 10;
                [vec_t, flag, RELRES, iter_num ] = gmres(mtx_A, rhs, RESTRT, linear_system_tol, pcg_iter_number, [], [], init_vec); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num(1)*RESTRT+iter_num(2), RELRES, flag);
            
            otherwise
                [vec_t, flag, RELRES, iter_num ] = bicgstabl(mtx_A, rhs, linear_system_tol, pcg_iter_number, [], [], init_vec); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        end
    end
    
else
    
    if isa(A,'function_handle')
        switch LS_solver
            case 'pcg'
%                 [vec_t, flag, RELRES, iter_num ] = pcg(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
%                      linear_system_tol, pcg_iter_number,
%                      @(x)Precond_M_LS_GEP(x, sigma, prec_M));
                     
                              [vec_t, flag, RELRES, iter_num ] = pcg(A, prec_M .\ rhs, ...
                     linear_system_tol, pcg_iter_number);
                 vec_t = prec_M .\ vec_t;
%                 [vec_t, flag, RELRES, iter_num ] = pcg(A, rhs, ...
%                      linear_system_tol, pcg_iter_number, prec_M, prec_M);
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
    
            case 'minres'
%                 [vec_t, flag, RELRES, iter_num ] = minres(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
%                      linear_system_tol, pcg_iter_number, @(x)Precond_M_LS_GEP(x, sigma, prec_M)); 
                [vec_t, flag, RELRES, iter_num ] = minres(fun_shift_Ar, rhs, ...
                     linear_system_tol, pcg_iter_number,prec_M, prec_M); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'bicgstab'
                [vec_t, flag, RELRES, iter_num ] = bicgstab(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
                     linear_system_tol, pcg_iter_number, @(x)Precond_M_LS_GEP(x, sigma, prec_M)); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'bicgstabl'
                [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
                     linear_system_tol, pcg_iter_number, @(x)Precond_M_LS_GEP(x, sigma, prec_M)); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'cgs'
                [vec_t, flag, RELRES, iter_num ] = cgs(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
                     linear_system_tol, pcg_iter_number, @(x)Precond_M_LS_GEP(x, sigma, prec_M)); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'gmres'
                RESTRT       = 10;
                [vec_t, flag, RELRES, iter_num ] = gmres(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
                     RESTRT, linear_system_tol, pcg_iter_number, @(x)Precond_M_LS_GEP(x, sigma, prec_M)); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num(1)*RESTRT+iter_num(2), RELRES, flag);
            
            otherwise
                [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)SEP_shift_mtx_prod_vec(vec, A, sigma ), rhs, ...
                     linear_system_tol, pcg_iter_number, @(x)Precond_M_LS_GEP(x, sigma, prec_M)); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        end
    else
        mtx_A = A - sigma * speye(size(A,1));
        switch LS_solver
            case 'pcg'
                [vec_t, flag, RELRES, iter_num ] = pcg(mtx_A, rhs, linear_system_tol, pcg_iter_number, @(x)Precond_M_LS_GEP(x, sigma, prec_M)); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
    
            case 'minres'
                [vec_t, flag, RELRES, iter_num ] = minres(mtx_A, rhs, linear_system_tol, pcg_iter_number, @(x)Precond_M_LS_GEP(x, sigma, prec_M)); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'bicgstab'
                [vec_t, flag, RELRES, iter_num ] = bicgstab(mtx_A, rhs, linear_system_tol, pcg_iter_number, @(x)Precond_M_LS_GEP(x, sigma, prec_M)); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'bicgstabl'
                [vec_t, flag, RELRES, iter_num ] = bicgstabl(mtx_A, rhs, linear_system_tol, pcg_iter_number, @(x)Precond_M_LS_GEP(x, sigma, prec_M)); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'cgs'
                [vec_t, flag, RELRES, iter_num ] = cgs(mtx_A, rhs, linear_system_tol, pcg_iter_number, @(x)Precond_M_LS_GEP(x, sigma, prec_M)); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
            case 'gmres'
                RESTRT       = 10;
                [vec_t, flag, RELRES, iter_num ] = gmres(mtx_A, rhs, RESTRT, linear_system_tol, pcg_iter_number, @(x)Precond_M_LS_GEP(x, sigma, prec_M)); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num(1)*RESTRT+iter_num(2), RELRES, flag);
            
            otherwise
                [vec_t, flag, RELRES, iter_num ] = bicgstabl(mtx_A, rhs, linear_system_tol, pcg_iter_number, @(x)Precond_M_LS_GEP(x, sigma, prec_M)); 
                fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        end
    end
   
end

%inner_iter_count = inner_iter_count + iter_num;
%outer_iter_count = outer_iter_count + 1;

end

% =========================================================================
%
function sol = Precond_M_LS_GEP(x, sigma, prec_M)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sol = prec_M(x, sigma);

end
