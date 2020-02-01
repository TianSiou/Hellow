function Call_Back_eigenmode(hobj,event,h)
    global ith_wave_vector
    global Pgraph
    global flag_eigenmode
    global hline
    rotate3d off

    c_mode = get( h.radio_mode(1), 'value' ); % 1: isotropic; 0: bi-isotropic;
    
    if isempty(flag_eigenmode) || flag_eigenmode == 0
        flag_eigenmode = 1;
    end
    if flag_eigenmode == 1
        hline = cb_em(h,'start');
        set(hobj,'string','run eigenmode');
        flag_eigenmode = -1
    else
        set(hline,'visible','off');
        set(hobj,'string','choose wave vector');
        Pgraph.ith_wave_vector=ith_wave_vector;
        if c_mode == 1
            gui_eigenmode2('start', Pgraph.eigenvec_array(1:end,:,ith_wave_vector), Pgraph.eigvalue_array(:,ith_wave_vector),Pgraph.sn);
        else
            gui_eigenmode2('start', Pgraph.eigenvec_array(1:end/2,:,ith_wave_vector), Pgraph.eigvalue_array(:,ith_wave_vector),Pgraph.sn);
        end
        flag_eigenmode = 1
    end
end