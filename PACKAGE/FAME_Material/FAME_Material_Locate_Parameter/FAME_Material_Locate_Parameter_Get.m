function Pmaterial = FAME_Material_Locate_Parameter_Get(file_name)
temp = 1;
eval(['Pmaterial=',file_name,'(temp,temp);']);
end