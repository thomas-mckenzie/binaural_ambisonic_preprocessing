function s = GetVoronoiPlotandSolidAng(azimuth, elevation, plotflag)

% Returns values of solid angles from arrays of azimuth and elevation
% angles in degrees.
%
% ========================
% Inputs:
%     azimuth:            array of azimuth angles in degrees
%
%     elevation           array of elevation angles in degrees
%
%     plotFlag:           Enter 1 to plot a voronoi sphere plot
% 
% Adapted from original code by Bruno Luong: https://uk.mathworks.com/matlabcentral/fileexchange/40989-voronoi-sphere
% ================ Tom McKenzie, University of York, 2018 ================
% Version 1.0 - 5/3/2018

[x, y, z] = sph2cart(deg2rad(azimuth), deg2rad(elevation), ones(size(azimuth)));

xyz = [x; y; z];
xyz = bsxfun(@rdivide, xyz, sqrt(sum(xyz.^2,1)));
n = length(azimuth);
[P, K, voronoiboundary] = voronoisphere(xyz);

% Compute Solid Angles
s = vcell_solidangle(P, K);
s = s / sum(s); % so that the sum of all the solid angles = 1.

%% Graphic
if plotflag == 1   
    
    f = figure;%'rend','painters','pos',[10 10 700 430]);
    clf(f);
    set(f,'Renderer','zbuffer');
    ax = axes('Parent', f);
    hold(ax, 'on');
    axis(ax,'equal');
    
    plot3(ax, xyz(1,:),xyz(2,:),xyz(3,:),'w.'); %plot white dots for the points
    
    for k = 1:n
        X = voronoiboundary{k};      
        
        for i = 1:length(X(:,1))
            for j = 1:length(X(1,:))
                if isnan(X(i,j))
                    X(i,j) = 0;
                end
            end
        end
        
        fill3(X(1,:),X(2,:),X(3,:),s(k),'Parent',ax,'EdgeColor','w');
    end
    
    c2 = colorbar;
    c2.Label.String = 'Solid Angle';% (x 10^{-2})';
    
    view([45 18])
    
    xlabel ( 'X axis' )
    ylabel ( 'Y axis' )
    zlabel ( 'Z axis' )
    
    grid on
    grid minor
    ax.MinorGridAlphaMode = 'manual';
    ax.MinorGridAlpha = 0.1;
    ax.FontSize = 14;
    
    colormap(flipud(winter));
    
    axis(ax,'equal');
    axis(ax,[-1 1 -1 1 -1 1]);
    
    xticks([-1 0 1]);
    yticks([-1 0 1]);
    zticks([-1 0 1]);
    
    set(gcf, 'Color', 'w');
    
end
end






function [Vertices, K, voronoiboundary] = voronoisphere(xyz, varargin)
% [Vertices, K, voronoiboundary] = voronoisphere(xyz)
%
% Compute the voronoi's diagram of points on the spheres S^2
%
% INPUT:
%   xyz is (3 x n) array, coordinates of n points in R^3
%   they all must be be normalized to 1 (i.e., belong to the 2-sphere)
%
% OUTPUTS:
%   - Vertices is (3 x m) array, coordinates of the vertices of voronoi diagram
%   - K is (n x 1) cell, each K{j} contains the indices of the voronoi cell
%       vertices correspond to xyz(:,j).
%       Vertices are counter-clockwise oriented when looking from outside.
%   - voronoiboundary is (n x 1) cell, each voronoiboundary{j} contains the
%     descretized spherical polygonal coordinates of the voronoi cell.
%     they are (3 x mj) arrays, mj are number of discretized points.
%
% voronoisphere(..., 'resolution', rrad) to provide the resolution of the
% boundary. RRAD is in radian. Default value is 0.0349 (~ 2 degrees).
%   
% Author: Bruno Luong <brunoluong@yahoo.com>
% Date creation: 28/March/2013
%                29/March/201": specify orientation
%
%   See also: voronoi, voronoin

%{
% Example
n = 100;
xyz = randn(3,n);
xyz = bsxfun(@rdivide, xyz, sqrt(sum(xyz.^2,1)));

[P, K, voronoiboundary] = voronoisphere(xyz);

%% Graphic
f = figure(1);
clf(f);
set(f,'Renderer','zbuffer');
ax = axes('Parent', f);
hold(ax, 'on');
set(ax, 'Color', 'w');
plot3(ax, xyz(1,:),xyz(2,:),xyz(3,:),'wo');
clmap = cool();
ncl = size(clmap,1);
for k = 1:n
    X = voronoiboundary{k};
    cl = clmap(mod(k,ncl)+1,:);
    fill3(X(1,:),X(2,:),X(3,:),cl,'Parent',ax,'EdgeColor','w');
end
axis(ax,'equal');
axis(ax,[-1 1 -1 1 -1 1]);
%}

% Get the resolution in radian
options = struct(varargin{:});
if isfield(options, 'resolution')
    resolution = options.resolution;
else
    resolution = 2*pi/180; % 2 degrees
end

%%
npnts = size(xyz,2);
T = convhull(xyz.');
nt = size(T,1);
E = [T(:,[1 2]); T(:,[2 3]); T(:,[3 1])]; % pair indexes
E = sort(E,2);
[~, ~, J] = unique(E, 'rows');
if ~all(accumarray(J,1)==2)
    error('Topology issue due to numerical precision')
end

% Which 2 seeds the edge vertices correspond? 
allids = repmat((1:nt).',[3 1]);
k = accumarray(J, (1:3*nt).', [], @(k) {k});
k = [k{:}];
vid = allids(k.');

% each row is 2 cell ids of the edge
cellofedge = E(k(1,:),:); % ne x 2
ne = size(cellofedge,1);
edges = repmat((1:ne).',[2 1]);
edgeofcell = accumarray(cellofedge(:),edges, [], @(e) {e});

% Center of the circumscribed Delaunay triangles
Vertices = Center(xyz, T);

% Build the geodesic arcs that link two vertices
nedges = size(vid,1);
edgearcs = cell(1,nedges);
for k=1:nedges
    edgearcs{k} = Arc(Vertices(:,vid(k,:)), resolution);
end

% Build the contour of the voronoi cells
vs = sort(vid,2);
voronoiboundary = cell(size(edgeofcell));
K = cell(size(edgeofcell));
for k = 1:npnts
    % ordering and orientation of the edges
    v = cycling_edge(edgeofcell{k}, vid);
    v = oriented_edge(v, Vertices, xyz(:,k));
    [~, loc] = ismember(sort(v,2), vs, 'rows');
    % joint the arcs
    X = edgearcs(loc);
    flip = v(:,1)~=vid(loc,1);
    X(flip) = cellfun(@fliplr, X(flip), 'unif', false);
    X = cellfun(@(x) x(:,1:end-1), X, 'unif', false); % remove duplicated points
    voronoiboundary{k} = cat(2, X{:});
    % Keep the indices of the voronoi's hull
    K{k} = v(:,1);
end

end % voronoisphere

%%
function G = Arc(AB, resolution)
% return an discretized arc between two points A and B
A = AB(:,1);
B = AB(:,2);
AxB = cross(A,B);
AdB = dot(A,B);
Ap = cross(AxB, A);
Ap = Ap/norm(Ap);
theta = atan2(sqrt(sum(AxB.^2,1)), AdB); % > 0
npnts = max(ceil(theta/resolution),2); % at least 2 points
theta = linspace(0, theta, npnts);
G = A*cos(theta) + Ap*sin(theta);
end % Arc

%%
function P = Center(xyz, T)
% Center of the circumscribed Delaunay triangles
XYZ = reshape(xyz(:,T),[3 size(T)]);
A = XYZ(:,:,1);
B = XYZ(:,:,2);
C = XYZ(:,:,3);
A = A-C;
B = B-C;
A2B = bsxfun(@times, sum(A.^2,1), B);
B2A = bsxfun(@times, sum(B.^2,1), A);
AxB = cross(A,B,1);
P = cross(A2B - B2A, AxB, 1);
P = C + bsxfun(@times,P,1./(2*sum(AxB.^2,1)));
nP = sqrt(sum(P.^2,1));
P = bsxfun(@times, P, 1./nP);
s = dot(AxB,C);
P = bsxfun(@times, P, sign(s));
end % Center

%%
function v = cycling_edge(edges, vertexes)
% Chain the edges in cycle
u = vertexes(edges,:).';
n = size(u, 2);
[~, ~, I] = unique(u);
I = reshape(I,[2 n]);
J = repmat(1:n,[2 1]);
if ~all(accumarray(I(:), 1) == 2)
    error('Topology issue due to numerical precision')
end
K = accumarray(I(:), J(:), [], @(x) {x});
K = [K{:}];
v = zeros([n 2]);
p = 0;
q = 1;
% chain the edges
for j = 1:n
    i = K(:,q);
    if i(1) == p
        p = i(2);
    else
        p = i(1);
    end    
    i = I(:,p);
    if i(1) == q
        v(j,:) = u([1 2],p);
        q = i(2);
    else
        v(j,:) = u([2 1],p);
        q = i(1);
    end
end % for-loop
end % cycling_edge

%%
function v = oriented_edge(v, P, xyz)
% Orient the edges counter-clockwise
Q = null(xyz.');
E = P(:,v([1:end 1],1));
xy = Q'*E;
a = (xy(1,1:end-1)-xy(1,2:end))*(xy(2,1:end-1)+xy(2,2:end))';
if xor(a < 0, det([xyz Q]) < 0) % Combine orientation and directness of Q
    v = rot90(v,2);
end
end % oriented_edge





function s = vcell_solidangle(P, K)
% s = vcell_solidangle(P, K)
%
% Compute the solid angles of All voronoi cell
%
% P is (3 x m) array of vertices, coordinates of the vertices of voronoi diagram
% K is (n x 1) cell, each K{j} contains the indices of the voronoi cell
%
% Restrictions:
% - P must be unit vectors
% - For each cell vertices must be counter-clockwise oriented when looking from outside
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Date creation: 16/July/2014
%
% See also: voronoisphere

% Turn it to false to maximize speed
check_unit = true;

if check_unit
    u = sum(P.^2, 1);
    if any(abs(u-1) > 1e-6)
        error('vcell_solidangle:Pnotnormalized', ...
            'vcell_solidangle: P must be unit vectors');
    end
end

s = cellfun(@(k) one_vcell_solidangle(P(:,k)), K);

end % vcell_solidangle

%%
function omega = one_vcell_solidangle(v)
% omega = one_vcell_solidangle(v)
% Compute the solid angle of spherical polygonal defined by v
% v is (3 x n) matrix, each column is the coordinates of the vertex (unit vector, not checked)
% "correctly" oriented
%
% Ref: A. van Oosterom, J. Strackee:? A solid angle of a plane triangle.? – IEEE Trans. Biomed. Eng. 30:2 (1983); 125–126. 
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Date creation: 16/July/2014
%
% See also: vcell_solidangle

n = size(v,2);
s = zeros(1,n);
for i = 2:n-1
    T = v(:,[1 i i+1]);
    num = det(T);
    denom = 1 + sum(sum(T .* T(:,[2 3 1]), 1), 2);
    s(i) = num / denom;
end

omega = atan(s);
omega = 2 * sum(omega);

end % one_vcell_solidangle

