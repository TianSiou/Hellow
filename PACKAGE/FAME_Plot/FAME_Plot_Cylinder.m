function h = FAME_Plot_Cylinder(a1,a2,a3, bot_center,top_center,radius)
    bot_center = bot_center*[a1';a2';a3'];
    top_center = top_center*[a1';a2';a3'];
    for i = 1:length(radius)
        [X_cy,Y_cy,Z_cy] = cylinder(radius(i),20);

        cylinder_len = norm(top_center(i,:) - bot_center(i,:));
        Z_cy(2,:) = cylinder_len;

        R = a2e3_rotate((top_center(i,:) - bot_center(i,:))/cylinder_len);

        theta   = linspace(0,2*pi,60);
        ph      = linspace(0,radius(i),2);
        [t,p]   = meshgrid(theta,ph);
        [X_top_disk,Y_top_disk,Z_top_disk] = pol2cart(t,p,0);
        Z_top_disk = Z_top_disk*ones(2,length(X_top_disk));
        X_bot_disk = X_top_disk;
        Y_bot_disk = Y_top_disk;
        Z_bot_disk = Z_top_disk + cylinder_len;
        
        X_top_disk_temp = [X_top_disk(1,:) X_top_disk(2,:)] ;
        Y_top_disk_temp = [Y_top_disk(1,:) Y_top_disk(2,:)] ;
        Z_top_disk_temp = [Z_top_disk(1,:) Z_top_disk(2,:)] ;
        X_bot_disk_temp = [X_bot_disk(1,:) X_bot_disk(2,:)] ;
        Y_bot_disk_temp = [Y_bot_disk(1,:) Y_bot_disk(2,:)] ;
        Z_bot_disk_temp = [Z_bot_disk(1,:) Z_bot_disk(2,:)] ;
        X_cy_temp = [X_cy(1,:) X_cy(2,:)] ;
        Y_cy_temp = [Y_cy(1,:) Y_cy(2,:)] ;
        Z_cy_temp = [Z_cy(1,:) Z_cy(2,:)] ;
        n_top_disk   = length(X_top_disk_temp);
        n_bot_disk   = length(X_bot_disk_temp);
        n_cy        = length(X_cy_temp);

        R_top_disk    = kron(R,eye(n_top_disk));
        R_bot_disk    = kron(R,eye(n_bot_disk));
        R_cy          = kron(R,eye(n_cy));
        temp_top_disk = R_top_disk \ [X_top_disk_temp';Y_top_disk_temp';Z_top_disk_temp'];
        X_top_disk    = [temp_top_disk(1:n_top_disk/2)';temp_top_disk(n_top_disk/2+1:n_top_disk)'] + bot_center(i,1);
        Y_top_disk    = [temp_top_disk(n_top_disk+1:n_top_disk+n_top_disk/2)';temp_top_disk(n_top_disk+n_top_disk/2+1:2*n_top_disk)'] + bot_center(i,2);
        Z_top_disk    = [temp_top_disk(2*n_top_disk+1:2*n_top_disk+n_top_disk/2)';temp_top_disk(2*n_top_disk+n_top_disk/2+1:3*n_top_disk)'] + bot_center(i,3);      
        temp_bot_disk = R_bot_disk \ [X_bot_disk_temp';Y_bot_disk_temp';Z_bot_disk_temp'];
        X_bot_disk    = [temp_bot_disk(1:n_bot_disk/2)';temp_bot_disk(n_bot_disk/2+1:n_bot_disk)'] + bot_center(i,1);
        Y_bot_disk    = [temp_bot_disk(n_bot_disk+1:n_bot_disk+n_bot_disk/2)';temp_bot_disk(n_bot_disk+n_bot_disk/2+1:2*n_bot_disk)'] + bot_center(i,2);
        Z_bot_disk    = [temp_bot_disk(2*n_bot_disk+1:2*n_bot_disk+n_bot_disk/2)';temp_bot_disk(2*n_bot_disk+n_bot_disk/2+1:3*n_bot_disk)'] + bot_center(i,3);     
        temp_cy       = R_cy \ [X_cy_temp';Y_cy_temp';Z_cy_temp'];
        X_cy          = [temp_cy(1:n_cy/2)';temp_cy(n_cy/2+1:n_cy)'] + bot_center(i,1);
        Y_cy          = [temp_cy(n_cy+1:n_cy+n_cy/2)';temp_cy(n_cy+n_cy/2+1:2*n_cy)'] + bot_center(i,2);
        Z_cy          = [temp_cy(2*n_cy+1:2*n_cy+n_cy/2)';temp_cy(2*n_cy+n_cy/2+1:3*n_cy)'] + bot_center(i,3);     

        h(i) = surf([X_cy,nan(2,1),X_top_disk,nan(2,1),X_bot_disk],[Y_cy,nan(2,1),Y_top_disk,nan(2,1),Y_bot_disk],[Z_cy,nan(2,1),Z_top_disk,nan(2,1),Z_bot_disk]);
    end
end

function R = a2e3_rotate(vec_x)  
    vec_x = reshape(vec_x,3,1);
    R = eye(3);
    for  i = 1:2
        a = vec_x(i);
        b = vec_x(3);
        
        [c,s,r] = givens_rotation(b, a);
        if r < 1e-8
            R_temp = eye(3,3);
            R_temp(i,i) = c;
            R_temp(i,3) = -s;
            R_temp(3,i) = s;
            R_temp(3,3) = c;
        else
            R_temp = eye(3,3);
            R_temp(i,i) = c/r;
            R_temp(i,3) = -s/r;
            R_temp(3,i) = s/r;
            R_temp(3,3) = c/r;
        end
        R = R_temp*R;
        vec_x = R_temp*vec_x;
    end
end

function [c,s,r] = givens_rotation(a, b)
    if b == 0
        c = sign(a);
        if(c == 0)
            c = 1.0; %Unlike other languages, MatLab's sign function returns 0 on input 0.
        end
        s = 0;
        r = abs(a);
    elseif a == 0
        c = 0;
        s = sign(b);
        r = abs(b);
    elseif abs(a) > abs(b)
        t = b/a;
        u = sign(a)*abs(sqrt(1+t*t));
        c = 1/u;
        s = c*t;
        r = a*u;
    else
        t = a/b;
        u = sign(b)*abs(sqrt(1+t*t));
        s = 1/u;
        c = s*t;
        r = b*u;
    end
end