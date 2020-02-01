function FAME_Material_DisplayInfo(material, lattice_constant)
    display('-----------Lattice Info.---------')
    fprintf('  Prototype : %s\n', material.Prototype)
    fprintf('  Lattice type : %s\n', material.lattice_type)
    fprintf('  Lattice constants : %s\n', lattice_constant)
    fprintf('  Size ratio : %s\n', material.Prototype)
end