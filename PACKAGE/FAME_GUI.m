function FAME_GUI
    clear  
    clc
    global flag_eigenmode
    flag_eigenmode = 0;
    
    FAME_Folder_Manager( )
    if ismac || isunix
        addpath PACKAGE/FAME_GUI/
        addpath PACKAGE/FAME_GUI/GUI_Call_Back_functions 
        addpath PACKAGE/FAME_GUI/GUI_eigenmode 
    else 
        addpath PACKAGE\FAME_GUI\
        addpath PACKAGE\FAME_GUI\GUI_Call_Back_functions 
        addpath PACKAGE\FAME_GUI\GUI_eigenmode 
    end
    
    FAME_GUI_Fig_Main % main figure script

    FAME_GUI_Fig_Info % information figure script

    FAME_GUI_Default  % default setting script
    
    FAME_GUI_Callback_Main % callback setting script for main figure