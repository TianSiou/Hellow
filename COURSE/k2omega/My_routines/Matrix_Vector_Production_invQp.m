function [ vec_y, relres, LS_iter, LS_time ] = Matrix_Vector_Production_invQp( vec_x, mode, mtx_palindromic, funhand, fun_invPp, option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function return the answer of linear system  %
%               y = Q(nu_0)^{-1} x                  %
% or                                                %
%               y = Q(nu_0)^{-T} x                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   global inner_it inner_relres inner_time inner_it_pre inner_relres_pre inner_time_pre
   
    if strcmp(mode,'transpose') || strcmp(mode,'trans') || strcmp(mode,'transp')
        LS_FUN = @(v) funhand.palindromic.Qpt(v,option.shift_palindromic);
    else
        LS_FUN = @(v) funhand.palindromic.Qp(v,option.shift_palindromic);
    end
   
    
    if option.isprecond
        PC_FUN = @(v) fun_invPp( v, mode );
        switch option.solver
            case 'gmres'
                ts = tic;
                [vec_y,flag_pre,relres_pre,LS_iter_pre] = feval(option.solver_pretreat,  LS_FUN, vec_x , 1e-10, option.maxit, PC_FUN,[], option.v0);
                LS_time_pre = toc(ts);
                ts = tic;
                [vec_y,flag,relres,LS_iter_gmres] = gmres( LS_FUN, vec_x , option.restartSize, option.tol, option.maxit, PC_FUN, [], vec_y);
                LS_time_gmres = toc(ts);
                LS_iter = (LS_iter_gmres(1)-1)*option.restartSize + LS_iter_gmres(2) + LS_iter_pre;
                LS_time = LS_time_pre + LS_time_gmres;
                
                inner_time = [inner_time,LS_time_gmres];
                inner_it_pre = [inner_it_pre,LS_iter_pre];
                inner_relres_pre = [inner_relres_pre,relres_pre];
                inner_time_pre = [inner_time_pre,LS_time_pre];
            otherwise
                ts = tic;
                [vec_y,flag,relres,LS_iter] = feval(option.solver,  LS_FUN, vec_x , option.tol, option.maxit, PC_FUN,[], option.v0);
                LS_time = toc(ts);
                inner_time = [inner_time,LS_time];
        end
    else
%         switch option.linear_system.solver
%             case 'gmres'
%                 [vec_y,flag,relres,LS_iter] = gmres( LS_FUN, vec_x , option.linear_system.restartSize, option.linear_system.tol, option.linear_system.maxit, [], [], option.linear_system.v0);
%             otherwise
%                 [vec_y,flag,relres,LS_iter] = feval(option.linear_system.solver,  LS_FUN, vec_x , option.linear_system.tol, option.linear_system.maxit, [], [], option.linear_system.v0);
%         end
    end
    
    inner_it = [inner_it,LS_iter];
    inner_relres = [inner_relres,relres];
    
    
    
    switch option.solver
        case 'gmres'
            if flag == 0
                fprintf('\n%s stopped at iteration %d to a solution with relative residual %e. (use time %f.)\n', option.solver_pretreat, LS_iter_pre, relres_pre, LS_time_pre);    
                fprintf('%s converged at iteration %d (inner iteration %d) to a solution with relative residual %e. (use time %f.)\n', option.solver, LS_iter_gmres(1), LS_iter_gmres(2), relres, LS_time);    
            else
                fprintf('\n%s stopped at iteration %d to a solution with relative residual %e. (use time %f.)\n', option.solver_pretreat, LS_iter_pre, relres_pre, LS_time_pre);    
                fprintf('%s stopped at iteration %d (inner iteration %d) to a solution with relative residual %e. (use time %f.)\n', option.solver, LS_iter(1), LS_iter(2), relres, LS_time);   
            end
        otherwise
            if flag == 0
                fprintf('\n%s converged at iteration %d to a solution with relative residual %e. (use time %f.)\n', option.solver, LS_iter, relres, LS_time);    
            else
                fprintf('\n%s stopped at iteration %d to a solution with relative residual %e. (use time %f.)\n', option.solver, LS_iter, relres, LS_time);   
            end
    end
end

% % function vec_y = Matrix_Vector_Production_invPp(vec_x, mode, grid_num, lattice_type, mtx_palindromic, Lambdas)
% function vec_y = Matrix_Vector_Production_invPp(vec_x, mode, mtx_palindromic, fun_Tp, fun_Tps)
%    vec_y = fun_Tps(vec_x);
%    if strcmp(mode,'transpose') || strcmp(mode,'trans') || strcmp(mode,'transp')
%        vec_y = mtx_palindromic.invDpT_mtx*vec_y;
%    else
%        vec_y = mtx_palindromic.invDp_mtx*vec_y;
%    end
%    vec_y = fun_Tp(vec_y);
% end


% function vec_y = Matrix_Vector_Production_Qp_hat(vec_x, nu_0, mode, freq, grid_num, lattice_type, mtx_palindromic, mtx_B, Lambdas)
% % This return the result of matrix-vector product for inv(Pp)Qp
% %     if nu_0 ~= 1
%         a_eps = sum(mtx_B.B_eps)/length(mtx_B.B_eps);
%         vec_y = (a_eps - mtx_B.B_eps).*vec_x;
%         coef  = 0.25*(freq^2)*(nu_0^2 - 2*nu_0 + 1);
%         vec_y = vec_x + coef*Matrix_Vector_Production_invPp(vec_y, mode, grid_num, lattice_type, mtx_palindromic, Lambdas);
% %     else
% %         vec_y = vec_x;
% %     end
% end
