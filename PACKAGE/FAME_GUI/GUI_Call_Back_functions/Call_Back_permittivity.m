function Call_Back_permittivity(hobj,event,h)
    
    [Popt] = FAME_GUI_Info_Get(h);
    
    permitt = str2num( get(hobj,'string') );
    
    switch Popt.material.material_type
        case {'isotropic' 'biisotropic'} % isotropic, biisotropic
            fool_check_isotropic(permitt,h)
        case 'anisotropic' % anisotropic
            fool_check_anisotropic(permitt,h)
    end
    
    FAME_GUI_Info_Set('permittivity', get(hobj,'string'), h);
end
function fool_check_isotropic(permitt,h)
    if size(permitt,1)~=1 ||  size(permitt,2)~=1
        warndlg('An error input! Please confirm the input permittivity as an scalar.','Warning!');
        set(h.edit_material(1),'string','13');
    end
end
function fool_check_anisotropic(permitt,h)
    if size(permitt,1)~=3 ||  size(permitt,2)~=3
        warndlg('An error input! Please confirm the input permittivity as a 3x3 tensor.','Warning!');
        set(h.edit_material(1),'string','[21.26 0 -14i; 0 16 0; 14i 0 21.26]');
    end
    
    temp= permitt - permitt';
    if norm( temp ) > 1e-12
        warndlg('The permittivity tensor must be a hermitian 3x3 tensor!','Warning!');
        set(h.edit_material(1),'string',mat2str((permitt+permitt')/2));
    end
end