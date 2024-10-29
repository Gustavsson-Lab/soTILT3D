function [transformedControlPoints] = GlobalQuadraticReverseMapping(controlPoints,tdata, weight_type)
%  Given a location in output space (x,y,z), find the weighted average
%  of the quadratic polynomials from all control points that influence (x,y,z) to
%  compute the corresponding location in input space (u,v,w).


nControlPoints = size(controlPoints,1);

if tdata.ndims_in == 2
    
    
    
elseif tdata.ndims_in == 3
    
    x = controlPoints(:,1);
    y = controlPoints(:,2);
    z = controlPoints(:,3);
    u = zeros(nControlPoints,1);
    v = zeros(nControlPoints,1);
    w = zeros(nControlPoints,1);
    
    for icp = 1:nControlPoints
        
        [u(icp),v(icp),w(icp)] = GetUV_GlobalQuadratic_3D(x(icp),y(icp),z(icp),tdata);
        
    end
    
    transformedControlPoints = [u,v,w];
    
end

end



