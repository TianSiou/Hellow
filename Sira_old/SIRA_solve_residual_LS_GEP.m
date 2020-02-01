function [ vec_t ] = SIRA_solve_residual_LS_GEP(A, B, sigma, rhs, LS_solver, pcg_iter_number, ...
    linear_system_tol, prec_M, init_vec)
%Solve the linear system
%        ( A - tau I ) vec_t = rhs.
%    

% global inner_iter_count outer_iter_count

if ( nargin <= 7 )
    switch LS_solver
        case 'pcg'
            [vec_t, flag, RELRES, iter_num ] = pcg(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
    
        case 'minres'
            [vec_t, flag, RELRES, iter_num ] = minres(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'bicgstab'
            [vec_t, flag, RELRES, iter_num ] = bicgstab(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'bicgstabl'
            [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'cgs'
            [vec_t, flag, RELRES, iter_num ] = cgs(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'gmres'
            RESTRT       = 10;
            [vec_t, flag, RELRES, iter_num ] = gmres(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                RESTRT, linear_system_tol, pcg_iter_number); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num(1)*RESTRT+iter_num(2), RELRES, flag);
            
        otherwise
            [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
    end
elseif ( nargin == 8 )
    switch LS_solver
        case 'pcg'
            [vec_t, flag, RELRES, iter_num ] = pcg(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, prec_M); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
    
        case 'minres'
            [vec_t, flag, RELRES, iter_num ] = minres(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, prec_M); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'bicgstab'
            [vec_t, flag, RELRES, iter_num ] = bicgstab(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, prec_M); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'bicgstabl'
            [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, prec_M); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'cgs'
            [vec_t, flag, RELRES, iter_num ] = cgs(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, prec_M); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'gmres'
            RESTRT       = 10;
            [vec_t, flag, RELRES, iter_num ] = gmres(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                RESTRT, linear_system_tol, pcg_iter_number, prec_M); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num(1)*RESTRT+iter_num(2), RELRES, flag);
        
        otherwise
            [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, prec_M); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
            
    end    
    
else
    switch LS_solver
        case 'pcg'
            [vec_t, flag, RELRES, iter_num ] = pcg(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, [], [], init_vec); 
            fprintf('iter = %4.0f; residual = %12.4e;  = %2.0f \n',iter_num, RELRES, flag);
    
        case 'minres'
            [vec_t, flag, RELRES, iter_num ] = minres(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, [], [], init_vec); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'bicgstab'
            [vec_t, flag, RELRES, iter_num ] = bicgstab(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, [], [], init_vec); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'bicgstabl'
            [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, [], [], init_vec); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'cgs'
            [vec_t, flag, RELRES, iter_num ] = cgs(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, [], [], init_vec); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
        case 'gmres'
            RESTRT       = 10;
            [vec_t, flag, RELRES, iter_num ] = gmres(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                RESTRT, linear_system_tol, pcg_iter_number, [], [], init_vec); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num(1)*RESTRT+iter_num(2), RELRES, flag);
            
        otherwise
            [vec_t, flag, RELRES, iter_num ] = bicgstabl(@(vec)GEP_shift_mtx_prod_vec(vec, A, B, sigma ), rhs, ...
                linear_system_tol, pcg_iter_number, [], [], init_vec); 
            fprintf('iter = %4.0f; residual = %12.4e; flag = %2.0f \n',iter_num, RELRES, flag);
        
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
