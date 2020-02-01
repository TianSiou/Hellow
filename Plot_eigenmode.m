mx = Par.mesh.edge_len(1);
my = Par.mesh.edge_len(2);
mz = Par.mesh.edge_len(3);
nx = Par.mesh.grid_num(1);
ny = Par.mesh.grid_num(2);
nz = Par.mesh.grid_num(3);
s = nx * ny * nz;
x =linspace(0, mx, nx);
y =linspace(0, my, ny);
z =linspace(0, mz, nz);
vec1 = Par.lattice.lattice_vec_a(:, 1);
vec2 = Par.lattice.lattice_vec_a(:, 2);
vec3 = Par.lattice.lattice_vec_a(:, 3);
start_pt = [0; 0; 0];
[X, Y, Z] = meshgrid(x, y, z);

xslice = [];
yslice = [];
zslice = [];
k = 2 * pi;

% range1 = 0.4205;
% range2 = 0.5758;
range1 = 0.4305;
range2 = 0.4701;

[col row] = find(omega_array > range1 * k & omega_array < range2 * k);
size(col, 1)
for idx = 1 : size(col, 1)
    figure(idx)
    [idx row(idx) col(idx) omega_array(col(idx),row(idx)) / 2 / pi]
    Eigenvector = ele_field{row(idx)}(:, col(idx));
    %         Vector = abs(Eigenvector(1 : s));
    %         Vector = abs(Eigenvector(s + 1 : 2 * s));
    %         Vector = abs(Eigenvector(2 * s + 1 : 3 * s));
    Vector = abs(Eigenvector(1 : s) + Eigenvector(s + 1 : 2 * s) + Eigenvector(2 * s + 1 : 3 * s));
    Vector = reshape(Vector, nx, ny, nz);
    W = zeros(ny, nx, nz);
    for t = 1 : nz
        W(:, :, t) = Vector(:, :, t)';
    end
    for t = linspace(0, mz, 30)
        yslice = t;
        h = slice( X, Y, Z, W, xslice, yslice, zslice, 'nearest');
        hold on
        h.FaceColor = 'interp';
        h.EdgeColor = 'none';
        alpha color
        alpha scaled
    end
    h.DiffuseStrength = 0.8;
%     colorbar('horiz')
    
    FAME_Plot_Parallelepiped(start_pt,vec1,vec2,vec3,Par.lattice.lattice_constant,'color', 'computational');
    title(['\fontsize{14}\omega = ', num2str(omega_array(col(idx),row(idx)) / 2 / pi, 5)]);
    hold off
    axis equal
    axis tight
    view([19 0])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'ztick',[])
    %         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.5, 0.2, 0.2]);
    %     saveas(gcf,[num2str(n_supercell(1)),num2str(ny),'-',num2str(idx)],'epsc');
end


function FAME_Plot_Parallelepiped(start_pt,vec1,vec2,vec3,lattice_constant,plot_mode, lattice_vec_mode)
vec1     = reshape(vec1,3,1);
vec2     = reshape(vec2,3,1);
vec3     = reshape(vec3,3,1);
start_pt = reshape(start_pt,3,1);
if strcmp(plot_mode,'color') == 1
    position_a(:,1) = reshape((start_pt+vec1),3,1);
    position_a(:,2) = reshape((start_pt+vec2),3,1);
    position_a(:,3) = reshape((start_pt+vec3),3,1);
    
    rate = 1.1;
    if isfield(lattice_constant,'Permutation') == 1 && strcmp(lattice_vec_mode,'computational') == 1
        %             plot3( [start_pt(1);position_a(1,1)], [start_pt(2);position_a(2,1)], [start_pt(3);position_a(3,1)], 'r-', 'LineWidth', 3);
        %             plot3( [start_pt(1);position_a(1,2)], [start_pt(2);position_a(2,2)], [start_pt(3);position_a(3,2)], 'g-', 'LineWidth', 3);
        %             plot3( [start_pt(1);position_a(1,3)], [start_pt(2);position_a(2,3)], [start_pt(3);position_a(3,3)], 'b-', 'LineWidth', 3);
        mArrow3([start_pt(1) start_pt(2) start_pt(3)],[position_a(1,1) position_a(2,1) position_a(3,1)], 'facealpha', 1, 'color', 'red', 'stemWidth', 0.02, 'tipWidth', 0.04);
        mArrow3([start_pt(1) start_pt(2) start_pt(3)],[position_a(1,2) position_a(2,2) position_a(3,2)], 'facealpha', 1, 'color', 'green', 'stemWidth', 0.02, 'tipWidth', 0.04);
        mArrow3([start_pt(1) start_pt(2) start_pt(3)],[position_a(1,3) position_a(2,3) position_a(3,3)], 'facealpha', 1, 'color', 'blue', 'stemWidth', 0.02, 'tipWidth', 0.04);
        
        text( rate*position_a(1,1), rate*position_a(2,1), rate*position_a(3,1), 'a_1', 'color', 'black', 'FontSize', 12);
        text( rate*position_a(1,2), rate*position_a(2,2), rate*position_a(3,2), 'a_2', 'color', 'black', 'FontSize', 12);
        text( rate*position_a(1,3), rate*position_a(2,3), rate*position_a(3,3), 'a_3', 'color', 'black', 'FontSize', 12);
    elseif strcmp(lattice_vec_mode,'original') == 1
        position_a_perm = position_a(:,lattice_constant.Permutation);
        
        %             plot3( [start_pt(1);position_a_perm(1,1)], [start_pt(2);position_a_perm(2,1)], [start_pt(3);position_a_perm(3,1)], 'r-', 'LineWidth', 2);
        %             plot3( [start_pt(1);position_a_perm(1,2)], [start_pt(2);position_a_perm(2,2)], [start_pt(3);position_a_perm(3,2)], 'g-', 'LineWidth', 2);
        %             plot3( [start_pt(1);position_a_perm(1,3)], [start_pt(2);position_a_perm(2,3)], [start_pt(3);position_a_perm(3,3)], 'b-', 'LineWidth', 2);
        mArrow3([start_pt(1) start_pt(2) start_pt(3)],[position_a_perm(1,1) position_a_perm(2,1) position_a_perm(3,1)], 'facealpha', 1, 'color', 'red', 'stemWidth', 0.02, 'tipWidth', 0.04);
        mArrow3([start_pt(1) start_pt(2) start_pt(3)],[position_a_perm(1,2) position_a_perm(2,2) position_a_perm(3,2)], 'facealpha', 1, 'color', 'green', 'stemWidth', 0.02, 'tipWidth', 0.04);
        mArrow3([start_pt(1) start_pt(2) start_pt(3)],[position_a_perm(1,3) position_a_perm(2,3) position_a_perm(3,3)], 'facealpha', 1, 'color', 'blue', 'stemWidth', 0.02, 'tipWidth', 0.04);
        
        %             text( rate*position_a(1,lattice_constant.Permutation(1)), rate*position_a(2,lattice_constant.Permutation(1)), rate*position_a(3,lattice_constant.Permutation(1)), '\bf{a}_1', 'color', 'black', 'FontSize', 14);
        %             text( rate*position_a(1,lattice_constant.Permutation(2)), rate*position_a(2,lattice_constant.Permutation(2)), rate*position_a(3,lattice_constant.Permutation(2)), '\bf{a}_2', 'color', 'black', 'FontSize', 14);
        %             text( rate*position_a(1,lattice_constant.Permutation(3)), rate*position_a(2,lattice_constant.Permutation(3)), rate*position_a(3,lattice_constant.Permutation(3)), '\bf{a}_3', 'color', 'black', 'FontSize', 14);
    end
elseif strcmp(plot_mode,'no_color') == 1
    plot3( [start_pt(1);start_pt(1)+vec1(1)'], [start_pt(2);start_pt(2)+vec1(2)'], [start_pt(3);start_pt(3)+vec1(3)'], 'k-',...
        [start_pt(1);start_pt(1)+vec2(1)'], [start_pt(2);start_pt(2)+vec2(2)'], [start_pt(3);start_pt(3)+vec2(3)'], 'k-',...
        [start_pt(1);start_pt(1)+vec3(1)'], [start_pt(2);start_pt(2)+vec3(2)'], [start_pt(3);start_pt(3)+vec3(3)'], 'k-');
end
vertex = [  start_pt(1) start_pt(2) start_pt(3);
    [start_pt(1) start_pt(2) start_pt(3)] + vec1';
    [start_pt(1) start_pt(2) start_pt(3)] + vec1'+vec2';
    [start_pt(1) start_pt(2) start_pt(3)] + vec2';
    [start_pt(1) start_pt(2) start_pt(3)] + vec3';
    [start_pt(1) start_pt(2) start_pt(3)] + vec3'+vec1';
    [start_pt(1) start_pt(2) start_pt(3)] + vec3'+vec1'+vec2';
    [start_pt(1) start_pt(2) start_pt(3)] + vec3'+vec2';
    [start_pt(1) start_pt(2) start_pt(3)] + 0.5*vec1';
    [start_pt(1) start_pt(2) start_pt(3)] + 0.5*vec1' + vec2';
    [start_pt(1) start_pt(2) start_pt(3)] + 0.5*vec1' + vec3';
    [start_pt(1) start_pt(2) start_pt(3)] + 0.5*vec1' + vec2' + vec3';];
edge = [1 2;
    2 3;
    3 4;
    4 1;
    1 5;
    2 6;
    3 7;
    4 8;
    5 6;
    6 7;
    7 8;
    8 5;
    9 10;
    9 11;
    12 10;
    12 11];
for i = 1:size(edge,1)
    idx1 = edge(i,1);
    idx2 = edge(i,2);
    path = [vertex(idx1,:);vertex(idx2,:)];
    if i == 1
        %             plot3( path(:,1), path(:,2), path(:,3), 'r-');
    elseif i == 4
        %             plot3( path(:,1), path(:,2), path(:,3), 'g-');
    elseif i == 5
        %             plot3( path(:,1), path(:,2), path(:,3), 'b-');
    else
        plot3( path(:,1), path(:,2), path(:,3), 'k-')
    end
    hold on
    %         if i == 1
    %             plot3( path(1,1), path(1,2), path(1,3), 'ro')
    %         end
end
end

function h = mArrow3(p1,p2,varargin)
%mArrow3 - plot a 3D arrow as patch object (cylinder+cone)
%
% syntax:   h = mArrow3(p1,p2)
%           h = mArrow3(p1,p2,'propertyName',propertyValue,...)
%
% with:     p1:         starting point
%           p2:         end point
%           properties: 'color':      color according to MATLAB specification
%                                     (see MATLAB help item 'ColorSpec')
%                       'stemWidth':  width of the line
%                       'tipWidth':   width of the cone                       
%
%           Additionally, you can specify any patch object properties. (For
%           example, you can make the arrow semitransparent by using
%           'facealpha'.)
%                       
% example1: h = mArrow3([0 0 0],[1 1 1])
%           (Draws an arrow from [0 0 0] to [1 1 1] with default properties.)
%
% example2: h = mArrow3([0 0 0],[1 1 1],'color','red','stemWidth',0.02,'facealpha',0.5)
%           (Draws a red semitransparent arrow with a stem width of 0.02 units.)
%
% hint:     use light to achieve 3D impression
%

propertyNames = {'edgeColor'};
propertyValues = {'none'};    

%% evaluate property specifications
for argno = 1:2:nargin-2
    switch varargin{argno}
        case 'color'
            propertyNames = {propertyNames{:},'facecolor'};
            propertyValues = {propertyValues{:},varargin{argno+1}};
        case 'stemWidth'
            if isreal(varargin{argno+1})
                stemWidth = varargin{argno+1};
            else
                warning('mArrow3:stemWidth','stemWidth must be a real number');
            end
        case 'tipWidth'
            if isreal(varargin{argno+1})
                tipWidth = varargin{argno+1};
            else
                warning('mArrow3:tipWidth','tipWidth must be a real number');
            end
        otherwise
            propertyNames = {propertyNames{:},varargin{argno}};
            propertyValues = {propertyValues{:},varargin{argno+1}};
    end
end            

%% default parameters
if ~exist('stemWidth','var')
    ax = axis;
    if numel(ax)==4
        stemWidth = norm(ax([2 4])-ax([1 3]))/300;
    elseif numel(ax)==6
        stemWidth = norm(ax([2 4 6])-ax([1 3 5]))/300;
    end
end
if ~exist('tipWidth','var')
    tipWidth = 3*stemWidth;
end
tipAngle = 22.5/180*pi;
tipLength = tipWidth/tan(tipAngle/2);
ppsc = 50;  % (points per small circle)
ppbc = 250; % (points per big circle)

%% ensure column vectors
p1 = p1(:);
p2 = p2(:);

%% basic lengths and vectors
x = (p2-p1)/norm(p2-p1); % (unit vector in arrow direction)
y = cross(x,[0;0;1]);    % (y and z are unit vectors orthogonal to arrow)
if norm(y)<0.1
    y = cross(x,[0;1;0]);
end
y = y/norm(y);
z = cross(x,y);
z = z/norm(z);

%% basic angles
theta = 0:2*pi/ppsc:2*pi; % (list of angles from 0 to 2*pi for small circle)
sintheta = sin(theta);
costheta = cos(theta);
upsilon = 0:2*pi/ppbc:2*pi; % (list of angles from 0 to 2*pi for big circle)
sinupsilon = sin(upsilon);
cosupsilon = cos(upsilon);

%% initialize face matrix
f = NaN([ppsc+ppbc+2 ppbc+1]);

%% normal arrow
if norm(p2-p1)>tipLength
    % vertices of the first stem circle
    for idx = 1:ppsc+1
        v(idx,:) = p1 + stemWidth*(sintheta(idx)*y + costheta(idx)*z);
    end
    % vertices of the second stem circle
    p3 = p2-tipLength*x;
    for idx = 1:ppsc+1
        v(ppsc+1+idx,:) = p3 + stemWidth*(sintheta(idx)*y + costheta(idx)*z);
    end
    % vertices of the tip circle
    for idx = 1:ppbc+1
        v(2*ppsc+2+idx,:) = p3 + tipWidth*(sinupsilon(idx)*y + cosupsilon(idx)*z);
    end
    % vertex of the tiptip
    v(2*ppsc+ppbc+4,:) = p2;

    % face of the stem circle
    f(1,1:ppsc+1) = 1:ppsc+1;
    % faces of the stem cylinder
    for idx = 1:ppsc
        f(1+idx,1:4) = [idx idx+1 ppsc+1+idx+1 ppsc+1+idx];
    end
    % face of the tip circle
    f(ppsc+2,:) = 2*ppsc+3:(2*ppsc+3)+ppbc;
    % faces of the tip cone
    for idx = 1:ppbc
        f(ppsc+2+idx,1:3) = [2*ppsc+2+idx 2*ppsc+2+idx+1 2*ppsc+ppbc+4];
    end

%% only cone v
else
    tipWidth = 2*sin(tipAngle/2)*norm(p2-p1);
    % vertices of the tip circle
    for idx = 1:ppbc+1
        v(idx,:) = p1 + tipWidth*(sinupsilon(idx)*y + cosupsilon(idx)*z);
    end
    % vertex of the tiptip
    v(ppbc+2,:) = p2;
    % face of the tip circle
    f(1,:) = 1:ppbc+1;
    % faces of the tip cone
    for idx = 1:ppbc
        f(1+idx,1:3) = [idx idx+1 ppbc+2];
    end
end

%% draw
fv.faces = f;
fv.vertices = v;
h = patch(fv);
for propno = 1:numel(propertyNames)
    try
        set(h,propertyNames{propno},propertyValues{propno});
    catch
        disp(lasterr)
    end
end
end