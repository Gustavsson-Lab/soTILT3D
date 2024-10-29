function [trans,uv,xy,uv_dev,xy_dev] = cp2tform_AG(varargin)
%CP2TFORM Infer spatial transformation from control point pairs.
%   CP2TFORM takes pairs of control points and uses them to infer a
%   spatial transformation.
%   This function was adapted from the MATLAB built-in function cp2tform().


%   TFORM = cp2tform_AG(cp_channel1,cp_channel2,tform_mode, nEqn, ...
%   weight_type, kthNeighbor, smoothnessParameter)
%   The local weighted mean (lwm) method creates a mapping, by inferring a polynomial at each
%   control point using neighboring control points. The mapping at any
%   location depends on a weighted average of these polynomials.  You can
%   optionally specify the number of points, nEqn, used to infer each
%   polynomial. The nEqn closest points are used to infer a polynomial
%   for each control point pair. nEqn can be as small as 10(?), BUT making N
%   small risks generating ill-conditioned polynomials.

%   tform_mode
%   -------------
%   cp2tform_AG requires a minimum number of control point pairs to infer a
%   TFORM structure of each tform_mode:
%
%       tform_mode            MINIMUM NUMBER OF PAIRS
%       -------------         -----------------------

%       'lwlinear'   (ORDER=1, 2D)       3
%       'lwlinear'   (ORDER=1, 3D)       4
%       'lwquadratic'(ORDER=2, 2D)       6
%       'lwquadratic'(ORDER=2, 3D)       10
%
%   When the minimum number of control points pairs are used for a particular
%   transformation, the coefficients are found exactly. If more than the minimum is used,
%   a least squares solution is found. See MLDIVIDE.
%
%   Note
%   ----
%   When either INPUT_POINTS or BASE_POINTS has a large offset with
%   respect to their origin (relative to range of values that it spans), the
%   points are shifted to center their bounding box on the origin before
%   fitting a TFORM structure.  This enhances numerical stability and
%   is handled transparently by wrapping the origin-centered TFORM within a
%   custom TFORM that automatically applies and undoes the coordinate shift
%   as needed. This means that fields(T) may give different results for
%   different coordinate inputs, even for the same TRANSFORMTYPE.


[uv, xy, method, weight_type, options] = ParseInputs(varargin{:});

% initialize deviation matrices
xy_dev = [];                %% transmitted channel, base points
uv_dev = [];                %% reflected channel, target points

% Assign function according to method and
% set K = number of control point pairs needed.
switch method
    case 'lwlinear'
        findT_fcn = @findLWLinear;
    case 'lwquadratic'
        findT_fcn = @findLWQuadratic;
    case 'globalquadratic'
        findT_fcn = @findGlobalQuadratic;
    otherwise
        error(message('images:cp2tform:internalProblem'))
end

% error if user enters too few control point pairs
M = size(uv,1);
if M<options.K
    error(message('images:cp2tform:rankError', options.K, method))
end

% get offsets to apply to before/after spatial transformation
uvShift = getShift(uv);
xyShift = getShift(xy);
needToShift = any([uvShift xyShift] ~= 0);

if ~needToShift
    % infer transform
    [trans, output] = findT_fcn(uv,xy,options);
else
    % infer transform for shifted data
    [tshifted, output] = findT_fcn(applyShift(uv,uvShift),...
        applyShift(xy,xyShift),options);
    
    % construct custom tform with tshifted between forward and inverse shifts
    tdata = struct('uvShift',uvShift,'xyShift',xyShift,'tshifted',tshifted);
    trans = maketform('custom',2,2,@fwd,@inverse,tdata);
end

trans.method = method;
trans.weight_type = weight_type;

function shift = getShift(points)
tol = 1e+3;
minPoints = min(points);
maxPoints = max(points);
center = (minPoints + maxPoints) / 2;
span = maxPoints - minPoints;
if (span(1) > 0 && abs(center(1))/span(1) > tol) ||...
        (span(2) > 0 && abs(center(2))/span(2) > tol)
    shift = center;
else
    shift = [0 0];
end

function shiftedPoints = applyShift(points,shift)
shiftedPoints = bsxfun(@minus, points, shift);

function points = undoShift(shiftedPoints,shift)
points = bsxfun(@plus, shiftedPoints, shift);

function x = fwd(u,t)
x = undoShift(tformfwd(applyShift(u,t.tdata.uvShift),t.tdata.tshifted),...
    t.tdata.xyShift);

function u = inverse(x,t)
u = undoShift(tforminv(applyShift(x,t.tdata.xyShift),t.tdata.tshifted),...
    t.tdata.uvShift);

function [trans,output] = findLWLinear(uv,xy,options)
% This function evaluates parameters for a 2D/3D locally weighted linear transformation
%
%In 2D:
% For a linear transformation:
% Lx(x,y) = Ax(x-xi) + Bx(y-yi) + ui
% Ly(x,y) = Ay(x-xi) + By(y-yi) + vi
%
% We need to find coefficient (A,B) such that Lx(x,y),Ly(x,y)
% evaluate to uk,vk at neigboring points k, respectively
%
% [ u1 v1 ]   [ (x1-xi) (y1-yi) 1 ]   [ Ax Ay ]
% [ u2 v2 ] = [ (x2-xi) (y2-xi) 1 ] * [ Bx By ]
% [  :  : ]   [    :       :    : ]   [ ui vi ]
% [ uk vk ]   [ (xk-xi) (yk-xi) 1 ]
%
%
%In 3D:
% For a linear transformation:
% Lx(x,y,z) = Ax(x-xi) + Bx(y-yi) + Cx(z-zi) + ui
% Ly(x,y,z) = Ay(x-xi) + By(y-yi) + Cy(z-zi) + vi
% Lz(x,y,z) = Az(x-xi) + Bz(y-yi) + Cz(z-zi) + wi
%
% We need to find coefficient (A,B,C) such that Lx(x,y,z),Ly(x,y,z),Lz(x,y,z)
% evaluate to uk,vk,wk at neigboring points k, respectively
%
% [ u1 v1 w1 ]   [ (x1-xi) (y1-yi) (z1-zi) 1 ]   [ Ax Ay Az ]
% [ u2 v2 w2 ] = [ (x2-xi) (y2-xi) (z2-zi) 1 ] * [ Bx By Ay ]
% [  :  :  : ]   [    :       :       :    : ]   [ Cx By Az ]
% [ uk vk wk ]   [ (xk-xi) (yk-xi) (zk-zi) 1 ]   [ ui vi wi ]
%
% Rewriting the above matrix equation:
% U = X * T
%
% With at least 4 correspondence points (3 neighbors), we can solve for T,
% T = X\U
%
%

output = [];
N = options.N;
N = options.nEqn;     % number of local correspondence points used to calculate the linear transformation
M = size(xy,1);     % number of correspondence points
k = options.kthNeighbor;              % specifies the sigma of the Gaussian weighting function as the kth closest point next to the control point
s = options.smoothnessParam;              % smoothness parameter - multiplies sigma(icp)

L = zeros(options.K,size(uv,2),M);
radii = zeros(M,1);
sigma = zeros(M,1);
% R = zeros(M,1);

if size(uv,2) == 2
    
    x = xy(:,1);
    y = xy(:,2);
    u = uv(:,1);
    v = uv(:,2);
    
    for icp = 1:M
        
        % find N closest points
        distcp = sqrt( (x-x(icp)).^2 + (y-y(icp)).^2 );
        [dist_sorted,indx] = sort(distcp);
        radii(icp) = dist_sorted(N);
        sigma(icp) = dist_sorted(k+1);                % the 1st closest control point
        % find the coordinates of the N closest points
        neighbors = indx(1:N);
        neighbors = sort(neighbors);
        xcp = x(neighbors);
        ycp = y(neighbors);
        ucp = u(neighbors);
        vcp = v(neighbors);
        
        % set up matrix eqn for local linear transformation
        % see Goshtasby, 2012, Eqn 9.97
        X = [xcp-x(icp),  ycp-y(icp), ones(N,1)];
        U = [ucp, vcp];
        
        if rank(X)>=options.K
            % options.K is the minimum number of points needed to solve the system of equations
            T = X\U;
            L(:,:,icp) = T;
        else
            error(message('images:cp2tform:rankError', options.K, 'polynomial'))
        end
        
    end
    
elseif size(uv,2) == 3
    
    x = xy(:,1);
    y = xy(:,2);
    z = xy(:,3);
    u = uv(:,1);
    v = uv(:,2);
    w = uv(:,3);
    
    for icp = 1:M
        
        % find N closest points
        distcp = sqrt( (x-x(icp)).^2 + (y-y(icp)).^2 + (z-z(icp)).^2);
        [dist_sorted,indx] = sort(distcp);
        %         radii(icp) = dist_sorted(N);
        sigma(icp) = dist_sorted(k+1);                % the 1st closest control point
        % find the coordinates of the N closest points
        neighbors = indx(1:N);
        neighbors = sort(neighbors);
        xcp = x(neighbors);
        ycp = y(neighbors);
        zcp = z(neighbors);
        ucp = u(neighbors);
        vcp = v(neighbors);
        wcp = w(neighbors);
        
        % set up matrix eqn for local linear transformation
        % see Goshtasby, 2012, Eqn 9.97
        X = [xcp-x(icp),  ycp-y(icp),  zcp-z(icp), ones(N,1)];
        U = [ucp, vcp, wcp];
        
        if rank(X)>=options.K
            % options.K is the minimum number of points needed to solve the system of equations
            T = X\U;
            L(:,:,icp) = T;
        else
            error(message('images:cp2tform:rankError', options.K, 'polynomial'))
        end
    end
    
end

tdata.LWLinearTData = L;
tdata.nEqn = N;
tdata.ControlPoints = xy;
tdata.RadiiOfInfluence = radii;
tdata.RadiusToKthNeighbor = sigma;
tdata.kthNeighbor = k;
tdata.smoothnessParam = s;

% Bybass the error checking done for the built-in matlab functions
% trans = maketform('custom',2,2,[],@inv_lwlinear,tdata);

trans.ndims_in = size(uv,2);
trans.ndims_out = size(uv,2);
trans.forward_fcn = [];
trans.inverse_fcn = @inv_lwlinear;
trans.tdata = tdata;

function [trans,output] = findLWQuadratic(uv,xy,options)
% This function evaluates parameters for a 2D/3D locally weighted quadratic transformation
%
% For a polynomial transformation:
%
% u = X*A, v = X*B, solve for A and B:
%     A = X\u;
%     B = X\v;
%
% The matrix X depends on the order of the polynomial.
% X will be M-by-K, where K = (order+1)*(order+2)/2;
% so A and B will be vectors of length K.
%
%   order = 2
%     X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2];
%     so X is an M-by-6 matrix
%
% see "Image registration by local approximation methods" Ardeshir
% Goshtasby, Image and Vision Computing, Vol 6, p. 255-261, 1988.


%In 2D:
% For a polynomial transformation:
% Lx(x,y) = Ax(x-xi) + Bx(y-yi) + Cx(x-xi)^2 + Dx(y-yi)^2 + Ex(x-xi)(y-yi) + ui
% Ly(x,y) = Ay(x-xi) + By(y-yi) + Cy(x-xi)^2 + Dy(y-yi)^2 + Ey(x-xi)(y-yi) + vi
%
% We need to find coefficient (A,B,C,D,E) such that Lx(x,y),Ly(x,y)
% evaluate to uk,vk at neigboring points k, respectively
%
% [ u1 v1 ]   [ (x1-xi) (y1-yi) (x1-xi)^2 (y1-yi)^2 (x1-xi)(y1-yi) 1 ]   [ Ax Ay ]
% [ u2 v2 ] = [ (x2-xi) (y2-xi) (x2-xi)^2 (y2-yi)^2 (x2-xi)(y2-yi) 1 ] * [ Bx By ]
% [  :  : ]   [    :       :       :          :            :       : ]   [ Cx Cy ]
% [ uk vk ]   [ (xk-xi) (yk-xi) (xk-xi)^2 (yk-yi)^2 (xk-xi)(yk-yi) 1 ]   [ Dx Dy ]
%                                                                        [ Ex Ey ]
%                                                                        [ ui vi ]
%

%In 3D:
% For a linear transformation:
% Lx(x,y,z) = Ax(x-xi) + Bx(y-yi) + Cx(z-zi) + Dx(x-xi)^2 + Ex(y-yi)^2 + Fx(z-zi)^2 + Gx(x-xi)(y-yi) + Hx(x-xi)(z-zi) + Ix(y-yi)(z-zi) + ui
% Ly(x,y,z) = Ay(x-xi) + By(y-yi) + Cy(z-zi) + Dy(x-xi)^2 + Ey(y-yi)^2 + Fy(z-zi)^2 + Gy(x-xi)(y-yi) + Hy(x-xi)(z-zi) + Iy(y-yi)(z-zi) + vi
% Lz(x,y,z) = Az(x-xi) + Bz(y-yi) + Cz(z-zi) + Dz(x-xi)^2 + Ez(y-yi)^2 + Fz(z-zi)^2 + Gz(x-xi)(y-yi) + Hz(x-xi)(z-zi) + Iy(y-yi)(z-zi) + wi
%
% We need to find coefficient (A,B,C,D,E,F,G,H,I) such that Lx(x,y,z),Ly(x,y,z),Lz(x,y,z)
% evaluate to uk,vk,wk at neigboring points k, respectively
%
% [ u1 v1 w1 ]   [ (x1-xi) (y1-yi) (z1-zi) (x1-xi)^2 (y1-yi)^2 (z1-zi)^2 (x1-xi)(y1-yi) (x1-xi)(z1-zi) (y1-yi)(z1-zi) 1 ]   [ Ax Ay Az ]
% [ u2 v2 w2 ] = [ (x2-xi) (y2-xi) (z2-zi) (x2-xi)^2 (y2-yi)^2 (z2-zi)^2 (x2-xi)(y2-yi) (x2-xi)(z2-zi) (y2-yi)(z2-zi) 1 ] * [ Bx By Ay ]
% [  :  :  : ]   [    :       :       :        :         :          :           :             :               :       : ]   [ Cx Bz Az ]
% [ uk vk wk ]   [ (xk-xi) (yk-xi) (zk-zi) (x1-xi)^2 (y1-yi)^2 (z1-zi)^2 (x1-xi)(y1-yi) (x1-xi)(z1-zi) (y1-yi)(z1-zi) 1 ]   [ Dx Dy Dz ]
%                                                                                                                           [ Ex Ey Ez ]
%                                                                                                                           [ Fx Fy Fz ]
%                                                                                                                           [ Gx Gy Gz ]
%                                                                                                                           [ Hx Hy Hz ]
%                                                                                                                           [ Ix Iy Iz ]
%                                                                                                                           [ ui vi wi ]
% Rewriting the above matrix equation:
% U = X * T
%
% With at least 10 correspondence points (9 neighbors), we can solve for T,
% T = X\U
%
% if (options.order ~= 2)
%     error(message('images:cp2tform:internalProblemPolyOrd'))
% end

output = [];
%N = options.N;
N = options.nEqn;     % number of local correspondence points used to calculate the linear transformation
M = size(xy,1);     % number of correspondence points
k = options.kthNeighbor;              % specifies the sigma of the Gaussian weighting function as the kth closest point next to the control point
s = options.smoothnessParam;              % smoothness parameter - multiplies sigma(icp)

L = zeros(options.K,size(uv,2),M);
radii = zeros(M,1);
sigma = zeros(M,1);
% R = zeros(M,1);

if size(uv,2) == 2
    
    x = xy(:,1);
    y = xy(:,2);
    u = uv(:,1);
    v = uv(:,2);
    
    for icp = 1:M
        
        % find N closest points
        distcp = sqrt( (x-x(icp)).^2 + (y-y(icp)).^2 );
        [dist_sorted,indx] = sort(distcp);
        radii(icp) = dist_sorted(N);
        sigma(icp) = dist_sorted(k+1);                % the 1st closest control point
        % find the coordinates of the N closest points
        neighbors = indx(1:N);
        neighbors = sort(neighbors);
        xcp = x(neighbors);
        ycp = y(neighbors);
        ucp = u(neighbors);
        vcp = v(neighbors);
        
        % set up matrix eqn for local linear transformation
        % see Goshtasby, 2012, Eqn 9.97
        X = [xcp-x(icp),  ycp-y(icp), (xcp-x(icp)).^2, (ycp-y(icp)).^2, (xcp-x(icp)).*(ycp-y(icp)), ones(N,1)];
        U = [ucp, vcp];
        
        if rank(X)>=options.K
            % options.K is the minimum number of points needed to solve the system of equations
            T = X\U;
            L(:,:,icp) = T;
        else
            error(message('images:cp2tform:rankError', options.K, 'polynomial'))
        end
        
    end
    
elseif size(uv,2) == 3
    
    x = xy(:,1);
    y = xy(:,2);
    z = xy(:,3);
    u = uv(:,1);
    v = uv(:,2);
    w = uv(:,3);
    
    for icp = 1:M
        
        % find N closest points
        distcp = sqrt( (x-x(icp)).^2 + (y-y(icp)).^2 + (z-z(icp)).^2);
        [dist_sorted,indx] = sort(distcp);
        %         radii(icp) = dist_sorted(N);
        radii(icp) = dist_sorted(1+1);                % the 1st closest control point
        sigma(icp) = dist_sorted(k+1);                % the kth closest control point
        % find the coordinates of the N closest points
        neighbors = indx(1:N);
        neighbors = sort(neighbors);
        xcp = x(neighbors);
        ycp = y(neighbors);
        zcp = z(neighbors);
        ucp = u(neighbors);
        vcp = v(neighbors);
        wcp = w(neighbors);
        
        % set up matrix eqn for local linear transformation
        % see Goshtasby, 2012, Eqn 9.97
        X = [xcp-x(icp),  ycp-y(icp),  zcp-z(icp), (xcp-x(icp)).^2, (ycp-y(icp)).^2, (zcp-z(icp)).^2, (xcp-x(icp)).*(ycp-y(icp)), (xcp-x(icp)).*(zcp-z(icp)), (ycp-y(icp)).*(zcp-z(icp)), ones(N,1)];
        U = [ucp, vcp, wcp];
        
        if rank(X)>=options.K
            % options.K is the minimum number of points needed to solve the system of equations
            T = X\U;
            L(:,:,icp) = T;
        else
            %T = X\U;
            %L(:,:,icp) = T;
            error(message('images:cp2tform:rankError', options.K, 'polynomial'))
        end
    end
    
end

tdata.LWQuadraticTData = L;
tdata.nEqn = N;
tdata.ControlPoints = xy;
tdata.RadiiOfInfluence = radii;
tdata.RadiusToKthNeighbor = sigma;
tdata.kthNeighbor = k;
tdata.smoothnessParam = s;

% Bybass the error checking done for the built-in matlab functions
% trans = maketform('custom',2,2,[],@inv_lwlinear,tdata);

trans.ndims_in = size(uv,2);
trans.ndims_out = size(uv,2);
trans.forward_fcn = [];
trans.inverse_fcn = @inv_lwquadratic;
trans.tdata = tdata;

function [trans,output] = findGlobalQuadratic(uv,xy,options)
% This function evaluates parameters for a 2D/3D global quadratic transformation
%
% For a polynomial transformation:
%
% u = X*A, v = X*B, solve for A and B:
%     A = X\u;
%     B = X\v;
%
% The matrix X depends on the order of the polynomial.
%
%   order = 2
%     X = [ones(M,1),  x,  y,  x.*y,  x.^2,  y.^2];
%     so X is an M-by-6 matrix
%
% see "Image registration by local approximation methods" Ardeshir
% Goshtasby, Image and Vision Computing, Vol 6, p. 255-261, 1988.


%In 2D:
% For a polynomial transformation:
% Lx(x,y) = Ax(x-xi) + Bx(y-yi) + Cx(x-xi)^2 + Dx(y-yi)^2 + Ex(x-xi)(y-yi) + ui
% Ly(x,y) = Ay(x-xi) + By(y-yi) + Cy(x-xi)^2 + Dy(y-yi)^2 + Ey(x-xi)(y-yi) + vi
%
% We need to find coefficient (A,B,C,D,E) such that Lx(x,y),Ly(x,y)
% evaluate to uk,vk at neigboring points k, respectively
%
% [ u1 v1 ]   [ (x1-xi) (y1-yi) (x1-xi)^2 (y1-yi)^2 (x1-xi)(y1-yi) 1 ]   [ Ax Ay ]
% [ u2 v2 ] = [ (x2-xi) (y2-xi) (x2-xi)^2 (y2-yi)^2 (x2-xi)(y2-yi) 1 ] * [ Bx By ]
% [  :  : ]   [    :       :       :          :            :       : ]   [ Cx Cy ]
% [ uk vk ]   [ (xk-xi) (yk-xi) (xk-xi)^2 (yk-yi)^2 (xk-xi)(yk-yi) 1 ]   [ Dx Dy ]
%                                                                        [ Ex Ey ]
%                                                                        [ ui vi ]
%

%In 3D:
% For a linear transformation:
% Lx(x,y,z) = Ax(x-xi) + Bx(y-yi) + Cx(z-zi) + Dx(x-xi)^2 + Ex(y-yi)^2 + Fx(z-zi)^2 + Gx(x-xi)(y-yi) + Hx(x-xi)(z-zi) + Ix(y-yi)(z-zi) + ui
% Ly(x,y,z) = Ay(x-xi) + By(y-yi) + Cy(z-zi) + Dy(x-xi)^2 + Ey(y-yi)^2 + Fy(z-zi)^2 + Gy(x-xi)(y-yi) + Hy(x-xi)(z-zi) + Iy(y-yi)(z-zi) + vi
% Lz(x,y,z) = Az(x-xi) + Bz(y-yi) + Cz(z-zi) + Dz(x-xi)^2 + Ez(y-yi)^2 + Fz(z-zi)^2 + Gz(x-xi)(y-yi) + Hz(x-xi)(z-zi) + Iy(y-yi)(z-zi) + wi
%
% We need to find coefficient (A,B,C,D,E,F,G,H,I) such that Lx(x,y,z),Ly(x,y,z),Lz(x,y,z)
% evaluate to uk,vk,wk at neigboring points k, respectively
%
% [ u1 v1 w1 ]   [ (x1-xi) (y1-yi) (z1-zi) (x1-xi)^2 (y1-yi)^2 (z1-zi)^2 (x1-xi)(y1-yi) (x1-xi)(z1-zi) (y1-yi)(z1-zi) 1 ]   [ Ax Ay Az ]
% [ u2 v2 w2 ] = [ (x2-xi) (y2-xi) (z2-zi) (x2-xi)^2 (y2-yi)^2 (z2-zi)^2 (x2-xi)(y2-yi) (x2-xi)(z2-zi) (y2-yi)(z2-zi) 1 ] * [ Bx By Ay ]
% [  :  :  : ]   [    :       :       :        :         :          :           :             :               :       : ]   [ Cx Bz Az ]
% [ uk vk wk ]   [ (xk-xi) (yk-xi) (zk-zi) (x1-xi)^2 (y1-yi)^2 (z1-zi)^2 (x1-xi)(y1-yi) (x1-xi)(z1-zi) (y1-yi)(z1-zi) 1 ]   [ Dx Dy Dz ]
%                                                                                                                           [ Ex Ey Ez ]
%                                                                                                                           [ Fx Fy Fz ]
%                                                                                                                           [ Gx Gy Gz ]
%                                                                                                                           [ Hx Hy Hz ]
%                                                                                                                           [ Ix Iy Iz ]
%                                                                                                                           [ ui vi wi ]
% Rewriting the above matrix equation:
% U = X * T
%
% With at least 10 correspondence points (9 neighbors), we can solve for T,
% T = X\U
%

output = [];
%N = options.N;

M = size(xy,1);       % number of correspondence points
N = M;                % number of local correspondence points used to calculate the linear transformation

L = zeros(options.K,size(uv,2));
% radii = zeros(M,1);
% sigma = zeros(M,1);
% R = zeros(M,1);

if size(uv,2) == 2
    
    x = xy(:,1);
    y = xy(:,2);
    u = uv(:,1);
    v = uv(:,2);
    
    %     for icp = 1:M
    %
    %         % find N closest points
    %         distcp = sqrt( (x-x(icp)).^2 + (y-y(icp)).^2 );
    %         [dist_sorted,indx] = sort(distcp);
    %         radii(icp) = dist_sorted(N);
    %         sigma(icp) = dist_sorted(k+1);                % the 1st closest control point
    %         % find the coordinates of the N closest points
    %         neighbors = indx(1:N);
    %         neighbors = sort(neighbors);
    %         xcp = x(neighbors);
    %         ycp = y(neighbors);
    %         ucp = u(neighbors);
    %         vcp = v(neighbors);
    %
    %     end
    
    % find the average position of the control points
    xMean = mean(x); yMean = mean(y);
    
    % set up matrix eqn for local linear transformation
    % see Goshtasby, 2012, Eqn 9.97
    X = [x-xMean,  y-yMean, (x-xMean).^2, (y-yMean).^2, (x-xMean).*(y-yMean), ones(N,1)];
    U = [u, v];
    
    if rank(X)>=options.K
        % options.K is the minimum number of points needed to solve the system of equations
        T = X\U;
        L = T;
    else
        error(message('images:cp2tform:rankError', options.K, 'polynomial'))
    end
    
    
elseif size(uv,2) == 3
    
    x = xy(:,1);
    y = xy(:,2);
    z = xy(:,3);
    u = uv(:,1);
    v = uv(:,2);
    w = uv(:,3);
    
    %     for icp = 1:M
    %
    %         % find N closest points
    %         distcp = sqrt( (x-x(icp)).^2 + (y-y(icp)).^2 + (z-z(icp)).^2);
    %         [dist_sorted,indx] = sort(distcp);
    %         %         radii(icp) = dist_sorted(N);
    %         sigma(icp) = dist_sorted(k+1);                % the 1st closest control point
    %         % find the coordinates of the N closest points
    %         neighbors = indx(1:N);
    %         neighbors = sort(neighbors);
    %         xcp = x(neighbors);
    %         ycp = y(neighbors);
    %         zcp = z(neighbors);
    %         ucp = u(neighbors);
    %         vcp = v(neighbors);
    %         wcp = w(neighbors);
    %
    %         % set up matrix eqn for local linear transformation
    %         % see Goshtasby, 2012, Eqn 9.97
    %         X = [xcp-x(icp),  ycp-y(icp),  zcp-z(icp), (xcp-x(icp)).^2, (ycp-y(icp)).^2, (zcp-z(icp)).^2, (xcp-x(icp)).*(ycp-y(icp)), (xcp-x(icp)).*(zcp-z(icp)), (ycp-y(icp)).*(zcp-z(icp)), ones(N,1)];
    %         U = [ucp, vcp, wcp];
    %
    %         if rank(X)>=options.K
    %             % options.K is the minimum number of points needed to solve the system of equations
    %             T = X\U;
    %             L(:,:,icp) = T;
    %         else
    %             error(message('images:cp2tform:rankError', options.K, 'polynomial'))
    %         end
    %     end
    
    % find the average position of the control points
    xMean = mean(x); yMean = mean(y);  zMean = mean(z);
    
    % set up matrix eqn for local linear transformation
    % see Goshtasby, 2012, Eqn 9.97
    X = [x-xMean,  y-yMean,  z-zMean, (x-xMean).^2, (y-yMean).^2, (z-zMean).^2, (x-xMean).*(y-yMean), (x-xMean).*(z-zMean), (y-yMean).*(z-zMean), ones(N,1)];
    U = [u, v, w];
    
    
    
    if rank(X)>=options.K
        % options.K is the minimum number of points needed to solve the system of equations
        T = X\U;
        L = T;
    else
        error(message('images:cp2tform:rankError', options.K, 'polynomial'))
    end
    
end

tdata.GlobalQuadraticTData = L;
% tdata.nEqn = N;
tdata.ControlPoints = xy;
% tdata.RadiiOfInfluence = radii;
% tdata.RadiusToKthNeighbor = sigma;
% tdata.kthNeighbor = k;
% tdata.smoothnessParam = s;

% Bybass the error checking done for the built-in matlab functions
% trans = maketform('custom',2,2,[],@inv_lwlinear,tdata);

trans.ndims_in = size(uv,2);
trans.ndims_out = size(uv,2);
trans.forward_fcn = [];
trans.inverse_fcn = @inv_globalquadratic;
trans.tdata = tdata;

function [uv, xy, method, weight_type, options] = ParseInputs(varargin)
% This function parses the inputs

% defaults
options.order = 3;
options.K = [];
N = [];

% iptchecknargin(5,7,nargin,mfilename);

% % figure out if syntax is
% % CP2TFORM(CPSTRUCT,TRANSFORMTYPE,...) or
% % CP2TFORM(INPUT_POINTS,BASE_POINTS,TRANSFORMTYPE,...)
%
% if isa(varargin{1},'struct')
%     % TRANS = CP2TFORM(CPSTRUCT,TRANSFORMTYPE)
%     % TRANS = CP2TFORM(CPSTRUCT,'polynomial',ORDER)
%     % TRANS = CP2TFORM(CPSTRUCT,'lwm',N)
%
%     iptchecknargin(5,6,nargin,mfilename);
%
%     [uv,xy] = cpstruct2pairs(varargin{1});
%     method = getMethod(varargin{2});
%
%     nargs_to_go = nargin - 2;
%     if nargs_to_go > 0
%         args = varargin(3:end);
%     end
%
% else
% TRANS = CP2TFORM(INPUT_POINTS,BASE_POINTS,TRANSFORMTYPE)
% TRANS = CP2TFORM(INPUT_POINTS,BASE_POINTS,'polynomial',ORDER)
% TRANS = CP2TFORM(INPUT_POINTS,BASE_POINTS,'lwm',N)

%iptchecknargin(6,7,nargin,mfilename);
narginchk(6,7)
uv = varargin{1};
xy = varargin{2};
method = getMethod(varargin{3});
nEqn = varargin{4};        % number of control Point included in the estimation of the linear fucntions
weight_type = varargin{5};          % type of weight functions used
kthNeighbor = varargin{6};
if weight_type == 'Gaussian'
    smoothnessParam = varargin{7};  % smoothness parameter to adjust width of all Gaussian weights
end


% end

if size(uv,2) < 2 || size(uv,2) > 3 || size(xy,2) < 2 || size(xy,2) > 3
    error(message('images:cp2tform:invalidControlPointMatrix'))
end

if size(uv,1) ~= size(xy,1)
    error(message('images:cp2tform:needSameNumControlPoints'))
end

switch method
    case 'lwlinear'
        order = 1;
        if size(uv,2) == 2
            options.K = (order+1)*(order+2)/2;
        elseif size(uv,2) == 3
            options.K = (order+1)*(order+2)/2 +1;
        end
        options.order = order;
        options.nEqn = nEqn;
        options.weight_type = weight_type;
        options.kthNeighbor = kthNeighbor;
        options.smoothnessParam = smoothnessParam;
        
        
        if isempty(N)
            % conservative default N protects user from ill-conditioned polynomials
            %             N = 2*options.K;
            N = nEqn;
            
        else
            % validate N
            if ~isnumeric(N) || numel(N)~=1 || rem(N,1)~=0 || N<options.K
                error(message('images:cp2tform:invalidInputN', options.K));
            end
        end
        options.N = N;
        
    case 'lwquadratic'
        order = 2;
        if size(uv,2) == 2
            options.K = 6;
        elseif size(uv,2) == 3
            options.K = 10;
        end
        options.order = order;
        options.nEqn = nEqn;
        options.weight_type = weight_type;
        options.kthNeighbor = kthNeighbor;
        options.smoothnessParam = smoothnessParam;
        
        
        if isempty(N)
            % conservative default N protects user from ill-conditioned polynomials
            %             N = 2*options.K;
            N = nEqn;
            
        else
            % validate N
            if ~isnumeric(N) || numel(N)~=1 || rem(N,1)~=0 || N<options.K
                error(message('images:cp2tform:invalidInputN', options.K));
            end
        end
        options.N = N;
        
    case 'globalquadratic'
        order = 2;
        if size(uv,2) == 2
            options.K = 6;
        elseif size(uv,2) == 3
            options.K = 10;
        end
        options.order = order;
        
        if isempty(N)
            % conservative default N protects user from ill-conditioned polynomials
            %             N = 2*options.K;
            N = size(uv,1);
            
        else
            % validate N
            if ~isnumeric(N) || numel(N)~=1 || rem(N,1)~=0 || N<options.K
                error(message('images:cp2tform:invalidInputN', options.K));
            end
        end
        options.N = N;
        
    otherwise
        error(message('images:cp2tform:internalProblem'))
        
end


%-------------------------------
% Function  getMethod
%
function method = getMethod(method_string)

method_string = lower(method_string);

% Figure out which method to use
methods = {'lwlinear', 'lwquadratic', 'globalquadratic'};
if ischar(method_string)
    indx = strmatch(method_string, methods);
    switch length(indx)
        case 0
            error(message('images:cp2tform:unrecognizedTransformType', method_string))
        case 1
            method = methods{indx};
        otherwise
            error(message('images:cp2tform:ambiguousTransformType', method_string))
    end
else
    error(message('images:cp2tform:transformTypeIsNotString'))
end

if strcmp(method,'linear conformal')
    method = 'nonreflective similarity';
end
