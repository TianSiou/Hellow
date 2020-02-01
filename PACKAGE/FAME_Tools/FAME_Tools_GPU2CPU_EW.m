function Freq_array = FAME_Tools_GPU2CPU_EW(file_name, graph_num, path_string, part_num)
    % Read file
    Freq_array = textread(file_name);
    % Reshape
    Freq_array = reshape(Freq_array,graph_num,length(Freq_array)/graph_num);
    % Unitized
    Freq_array = sqrt(Freq_array);
    % Plot band structure
    FAME_Plot_Band_Structure( path_string, part_num, Freq_array )
end