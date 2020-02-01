function Call_Back_permeability(hobj,event,h)

    [Popt] = FAME_GUI_Info_Get(h);
    
    permeab = str2num( get(hobj,'string') );
    
    switch Popt.material.material_type
        case {'isotropic' 'biisotropic'} % isotropic, biisotropic
            fool_check_isotropic(permeab,h)
        case 'anisotropic' % anisotropic
            fool_check_anisotropic(permeab,h)
    end
    
    FAME_GUI_Info_Set('permeability', get(hobj,'string'), h);
end
function fool_check_isotropic(permeab,h)
    if size(permeab,1)~=1 ||  size(permeab,2)~=1
        warndlg('An error input! Please confirm the input permeability as an scalar.','Warning!');
        set(h.edit_material(2),'string','1');
    end
end
function fool_check_anisotropic(permeab,h)
    if size(permeab,1)~=3 ||  size(permeab,2)~=3
        warndlg('An error input! Please confirm the input permeability as a 3x3 tensor.','Warning!');
        set(h.edit_material(2),'string','[1 0 0; 0 1 0; 0 0 1]');
    end
end