function Fun_InverseWoodpile = FAME_Material_Locate_Handle_ORC_InverseWoodpile( lattice_constant, cylinder_radius )
    
    temp_mod_c = @(x) mod( x, lattice_constant.c );
    temp_mod_a = @(z) mod( z, lattice_constant.a );

%     r_1 = @(x,y,z) reshape( sqrt( (temp_mod_a(x) -                      cylinder_radius).^2 + (temp_mod_c(z)                         ).^2 ), 1, length(x));
%     r_2 = @(x,y,z) reshape( sqrt( (temp_mod_a(x) -                      cylinder_radius).^2 + (temp_mod_c(z) -    lattice_constant.c ).^2 ), 1, length(x));
%     r_3 = @(x,y,z) reshape( sqrt( (temp_mod_a(x) - lattice_constant.a + cylinder_radius).^2 + (temp_mod_c(z) - .5*lattice_constant.c ).^2 ), 1, length(x));
    
    r_1 = @(x,y,z) reshape( sqrt( (temp_mod_a(x) - 0.25*lattice_constant.a).^2 + (temp_mod_c(z)                         ).^2 ), 1, length(x));
    r_2 = @(x,y,z) reshape( sqrt( (temp_mod_a(x) - 0.25*lattice_constant.a).^2 + (temp_mod_c(z) -    lattice_constant.c ).^2 ), 1, length(x));
    r_3 = @(x,y,z) reshape( sqrt( (temp_mod_a(x) - 0.75*lattice_constant.a).^2 + (temp_mod_c(z) - .5*lattice_constant.c ).^2 ), 1, length(x));
    
    r_4 = @(x,y,z) reshape( sqrt( (temp_mod_c(y) - 0.5*lattice_constant.c).^2 + (temp_mod_a(x)                         ).^2 ), 1, length(x));
    r_5 = @(x,y,z) reshape( sqrt( (temp_mod_c(y) - 0.5*lattice_constant.c).^2 + (temp_mod_a(x) -     lattice_constant.a).^2 ), 1, length(x));
    r_6 = @(x,y,z) reshape( sqrt( (temp_mod_c(y)                         ).^2 + (temp_mod_a(x) - 0.5*lattice_constant.a).^2 ), 1, length(x));
    r_7 = @(x,y,z) reshape( sqrt( (temp_mod_c(y) -     lattice_constant.c).^2 + (temp_mod_a(x) - 0.5*lattice_constant.a).^2 ), 1, length(x));
    
    R = @(x,y,z) min( [ r_1(x,y,z); r_2(x,y,z); r_3(x,y,z); r_4(x,y,z); r_5(x,y,z); r_6(x,y,z); r_7(x,y,z)] );
    
    Fun_InverseWoodpile{1} = @(x,y,z,a1,a2,a3) R(x,y,z);
end