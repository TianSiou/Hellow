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
    
    solve_LS_PQEP = @(vec_x, transp_flag) Matrix_Vector_Production_invQp( vec_x, transp_flag, mtx.mtx_palindromic, funhand, fun_invPp, option.linear_system);
 
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