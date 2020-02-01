function idx = FAME_Material_Locate_Handle_No191_Hexagonal_temp( X, Y, Z, a1, a2, a3, ring_thinkness, ring_inner_cylinder_radius, ring_outer_cylinder_radius, connect_cylinder_radius )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See detail on following website:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	ring_cylinder_top_centers = [ 0, 0, ring_thinkness;
				                  0, 0, 1             ;
                                  1, 0, ring_thinkness;
				                  1, 0, 1             ;
                                  0, 1, ring_thinkness;
				                  0, 1, 1             ;
                                  1, 1, ring_thinkness;
				                  1, 1, 1             ]*[a1,a2,a3]';
	ring_cylinder_bot_centers = [ 0, 0, 0               ;
				                  0, 0, 1-ring_thinkness;
                                  1, 0, 0               ;
				                  1, 0, 1-ring_thinkness;
                                  0, 1, 0               ;
				                  0, 1, 1-ring_thinkness;
                                  1, 1, 0               ;
				                  1, 1, 1-ring_thinkness;]*[a1,a2,a3]';
	
    Point_idx_ring = [];
    for i = 1:size(ring_cylinder_top_centers,1)/2
        [ Point_idx_ring_inner_upper ] = FAME_Material_Locate_Cylinder( [X,Y,Z], ring_inner_cylinder_radius, ring_cylinder_bot_centers(2*i,:), ring_cylinder_top_centers(2*i,:) );
        [ Point_idx_ring_inner_lower ] = FAME_Material_Locate_Cylinder( [X,Y,Z], ring_inner_cylinder_radius, ring_cylinder_bot_centers(2*i-1,:), ring_cylinder_top_centers(2*i-1,:) );
        [ Point_idx_ring_outer_upper ] = FAME_Material_Locate_Cylinder( [X,Y,Z], ring_outer_cylinder_radius, ring_cylinder_bot_centers(2*i,:), ring_cylinder_top_centers(2*i,:) );
        [ Point_idx_ring_outer_lower ] = FAME_Material_Locate_Cylinder( [X,Y,Z], ring_outer_cylinder_radius, ring_cylinder_bot_centers(2*i-1,:), ring_cylinder_top_centers(2*i-1,:) );
        Point_idx_ring = union( Point_idx_ring, setdiff( union(Point_idx_ring_outer_upper,Point_idx_ring_outer_lower), union(Point_idx_ring_inner_upper,Point_idx_ring_inner_lower)) );
    end

	conncet_cylinder_top_centers = [ .5*(ring_outer_cylinder_radius + ring_inner_cylinder_radius), 0, 1;
                                                           1-.5*(ring_outer_cylinder_radius + ring_inner_cylinder_radius), 0, 1 ];
	conncet_cylinder_bot_centers = [ .5*(ring_outer_cylinder_radius + ring_inner_cylinder_radius), 0, 0;
                                                           1-.5*(ring_outer_cylinder_radius + ring_inner_cylinder_radius), 0, 0  ];
	conncet_cylinder_top_centers = [conncet_cylinder_top_centers; conncet_cylinder_top_centers + [0,1,0]];
	conncet_cylinder_bot_centers = [conncet_cylinder_bot_centers; conncet_cylinder_bot_centers + [0,1,0]];
	conncet_cylinder_top_centers = [conncet_cylinder_top_centers; conncet_cylinder_top_centers(:,[2,1,3])];
	conncet_cylinder_bot_centers = [conncet_cylinder_bot_centers; conncet_cylinder_bot_centers(:,[2,1,3])];
    conncet_cylinder_top_centers = [conncet_cylinder_top_centers; 0.45, 0.45, 0; 0.55, 0.55, 0];
    conncet_cylinder_bot_centers = [conncet_cylinder_bot_centers; 0.45, 0.45, 1; 0.55, 0.55, 1];
Point_idx_connect = [];
for i = 1:size(conncet_cylinder_top_centers,1)
    top_center = conncet_cylinder_top_centers(i,:)*[a1,a2,a3]';
    bot_center = conncet_cylinder_bot_centers(i,:)*[a1,a2,a3]';
	idx_temp = FAME_Material_Locate_Cylinder( [X,Y,Z], connect_cylinder_radius, top_center, bot_center );
	Point_idx_connect = union(Point_idx_connect, idx_temp);
end

idx = union(Point_idx_ring,Point_idx_connect);


		
end