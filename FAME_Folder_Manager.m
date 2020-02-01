function FAME_Folder_Manager( fame_m_path_string )
if exist(fame_m_path_string) ~= 7
    fame_m_path_string = pwd;
end
% addpath('C:\Users\User\Desktop\FAME_m_g071708\JD_SEP');
% addpath('C:\Users\User\Desktop\FAME_m_g071708\JD_SEP\precondition');
% addpath('C:\Users\User\Desktop\FAME_m_g071708\JD_SEP\setting');
% addpath('C:\Users\User\Desktop\FAME_m_g071708\JD_SEP\update');
addpath('C:\Users\User\Desktop\FAME_m_g071708\SIRA');
%% First layer folder
first_layer_folder_name = {'PACKAGE'};
for i = 1:length(first_layer_folder_name)
    path_string = fullfile( fame_m_path_string, first_layer_folder_name{i});
    addpath( path_string );
    fprintf('Add path: %s\n',path_string);
end
%% Second layer folder in PACKAGE
second_layer_folder_name = {'FAME_Tools',...
                            'FAME_Matrix',...
                            'FAME_Parameter',...
                            'FAME_Material',...
                            'FAME_Matrix_Vector_Production',...
                            'FAME_Fast_Algorithms',...
                            'FAME_Plot',...
                            'FAME_GUI',...
                            'FAME_fem'};
for i = 1:length(second_layer_folder_name)
    path_string = fullfile( fame_m_path_string, 'PACKAGE',second_layer_folder_name{i});
    addpath( path_string );
    fprintf('Add path: %s\n',path_string);
end
%% Third layer folder
third_layer_folder_name_fem = {'FAME_Anisotropic_FEM',...
                               'FAME_Biisotropic_FEM'};
third_layer_folder_name_GUI = {'GUI_Call_Back_functions',...
                               'GUI_Call_Back_functions_user_defined',...
                               'GUI_eigenmode',...
                               'GUI_User_Defined'};
third_layer_folder_name_Material = {'FAME_Material_Isofunction',...
                                    'FAME_Material_Locate_Handle',...
                                    'FAME_Material_Locate_Parameter'};
third_layer_folder_name_Matrix_Vector_Production = {'FAME_Matrix_Vector_Production_Isotropic',...
                                                    'FAME_Matrix_Vector_Production_Anisotropic',...
                                                    'FAME_Matrix_Vector_Production_Biisotropic',...
                                                    'FAME_Matrix_Vector_Production_FFT',...
                                                    'FAME_Matrix_Vector_Production_IFFT',...
                                                    'FAME_Matrix_Vector_Production_PML'};
% third_layer_folder_name_Fast_Algorithms = {};
% third_layer_folder_name_Plot = {};
% third_layer_folder_name_Tools = {};
% third_layer_folder_name_Matrix = {};
% third_layer_folder_name_Parameter = {};


for i = 1:length(third_layer_folder_name_fem)
    path_string = fullfile( fame_m_path_string, 'PACKAGE','FAME_fem',third_layer_folder_name_fem{i});
    addpath( path_string );
    fprintf('Add path: %s\n',path_string);
end                                                               
for i = 1:length(third_layer_folder_name_GUI)
    path_string = fullfile( fame_m_path_string, 'PACKAGE','FAME_GUI',third_layer_folder_name_GUI{i});
    addpath( path_string );
    fprintf('Add path: %s\n',path_string);
end                                                        
for i = 1:length(third_layer_folder_name_Material)
    path_string = fullfile( fame_m_path_string, 'PACKAGE','FAME_Material',third_layer_folder_name_Material{i});
    addpath( path_string );
    fprintf('Add path: %s\n',path_string);
end                                                       
for i = 1:length(third_layer_folder_name_Matrix_Vector_Production)
    path_string = fullfile( fame_m_path_string, 'PACKAGE','FAME_Matrix_Vector_Production',third_layer_folder_name_Matrix_Vector_Production{i});
    addpath( path_string );
    fprintf('Add path: %s\n',path_string);
end                  



% if ismac || isunix
% %% MAC/Linus system   
%     addpath PACKAGE
%     addpath PACKAGE/FAME_Tools
%     addpath PACKAGE/FAME_Matrix
%     addpath PACKAGE/FAME_Parameter
%     addpath PACKAGE/FAME_Material
%     addpath PACKAGE/FAME_Material/FAME_Material_Locate_Handle
%     addpath PACKAGE/FAME_Material/FAME_Material_Locate_Parameter
%     addpath PACKAGE/FAME_Material/FAME_Material_Isofunction
%     addpath PACKAGE/FAME_Matrix_Vector_Production
%     addpath PACKAGE/FAME_Fast_Algorithms
%     addpath PACKAGE/FAME_Plot
%     addpath PACKAGE/FAME_Material
%     addpath PACKAGE/FAME_Matrix_Vector_Production/FAME_Matrix_Vector_Production_Isotropic
%     addpath PACKAGE/FAME_Matrix_Vector_Production/FAME_Matrix_Vector_Production_Anisotropic
%     addpath PACKAGE/FAME_Matrix_Vector_Production/FAME_Matrix_Vector_Production_Biisotropic
%     addpath PACKAGE/FAME_Matrix_Vector_Production/FAME_Matrix_Vector_Production_PML
%     addpath PACKAGE/FAME_Matrix_Vector_Production/FAME_Matrix_Vector_Production_FFT
%     addpath PACKAGE/FAME_Matrix_Vector_Production/FAME_Matrix_Vector_Production_IFFT
%     
%     addpath PACKAGE/FAME_GUI/
%     addpath PACKAGE/FAME_GUI/GUI_eigenmode 
% 
% else
% %% Windows system
%     addpath PACKAGE
%     addpath PACKAGE\FAME_Tools
%     addpath PACKAGE\FAME_Matrix
%     addpath PACKAGE\FAME_Parameter
%     addpath PACKAGE\FAME_Material
%     addpath PACKAGE\FAME_Material\FAME_Material_Locate_Handle
%     addpath PACKAGE\FAME_Material\FAME_Material_Locate_Parameter
%     addpath PACKAGE\FAME_Material\FAME_Material_Isofunction
%     addpath PACKAGE\FAME_Matrix_Vector_Production
%     addpath PACKAGE\FAME_Fast_Algorithms
%     addpath PACKAGE\FAME_Plot
%     addpath PACKAGE\FAME_Material
%     addpath PACKAGE\FAME_Matrix_Vector_Production\FAME_Matrix_Vector_Production_Isotropic
%     addpath PACKAGE\FAME_Matrix_Vector_Production\FAME_Matrix_Vector_Production_Anisotropic
%     addpath PACKAGE\FAME_Matrix_Vector_Production\FAME_Matrix_Vector_Production_Biisotropic
%     addpath PACKAGE\FAME_Matrix_Vector_Production\FAME_Matrix_Vector_Production_PML
%     addpath PACKAGE\FAME_Matrix_Vector_Production\FAME_Matrix_Vector_Production_FFT
%     addpath PACKAGE\FAME_Matrix_Vector_Production\FAME_Matrix_Vector_Production_IFFT
%     addpath PACKAGE\FAME_fem\FAME_Anisotropic_FEM
%     addpath PACKAGE\FAME_fem\FAME_Biisotropic_FEM
%     addpath PACKAGE\FAME_fem\
%     
%     addpath PACKAGE\FAME_GUI\
%     addpath PACKAGE\FAME_GUI\GUI_eigenmode 
% end
end