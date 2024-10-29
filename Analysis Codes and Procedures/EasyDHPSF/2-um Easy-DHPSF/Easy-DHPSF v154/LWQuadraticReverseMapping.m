function [transformedControlPoints] = LWQuadraticReverseMapping(controlPoints,tdata, weight_type)
%  Given a location in output space (x,y,z), find the weighted average
%  of the quadratic polynomials from all control points that influence (x,y,z) to
%  compute the corresponding location in input space (u,v,w).

nControlPoints = size(controlPoints,1);

if tdata.ndims_in == 2
    
    x = controlPoints(:,1);
    y = controlPoints(:,2);
    u = zeros(nControlPoints,1);
    v = zeros(nControlPoints,1);
    
    for icp = 1:nControlPoints
        
        Weights = zeros(nControlPoints,1);
        U = zeros(nControlPoints,1);
        V = zeros(nControlPoints,1);
        
        u_numerator = 0.0;
        v_numerator = 0.0;
        denominator = 0.0;
        
        [Weights, coordinates] = GetWeights_2D(x(icp),y(icp),tdata,weight_type);
        [U,V] = GetUV_LWQuadratic_2D(x(icp),y(icp),tdata,Weights);
        coordinates = find(~isnan(U));
        
        u_numerator = sum(Weights(coordinates).*U(coordinates));
        v_numerator = sum(Weights(coordinates).*V(coordinates));
        denominator = sum(Weights(coordinates));
        
        if denominator ~= 0.0
            u(icp) = u_numerator/denominator;
            v(icp) = v_numerator/denominator;
        else
            u(icp) = NaN;
            v(icp) = NaN;
            warning('no control points influence this (x,y)')
        end
    end
    clear W U V
    transformedControlPoints = [u,v];
    
elseif tdata.ndims_in == 3
    
    x = controlPoints(:,1);
    y = controlPoints(:,2);
    z = controlPoints(:,3);
    u = zeros(nControlPoints,1);
    v = zeros(nControlPoints,1);
    w = zeros(nControlPoints,1);
    
    for icp = 1:nControlPoints
        
        Weights = zeros(nControlPoints,1);
        U = zeros(nControlPoints,1);
        V = zeros(nControlPoints,1);
        W = zeros(nControlPoints,1);        
        
        u_numerator = 0.0;
        v_numerator = 0.0;
        w_numerator = 0.0;
        denominator = 0.0;
                
        [Weights, coordinates] = GetWeights_3D(x(icp),y(icp),z(icp),tdata,weight_type);
        [U,V,W] = GetUV_LWQuadratic_3D(x(icp),y(icp),z(icp),tdata,Weights);
        coordinates = find(~isnan(U));
        
        u_numerator = sum(Weights(coordinates).*U(coordinates));
        v_numerator = sum(Weights(coordinates).*V(coordinates));
        w_numerator = sum(Weights(coordinates).*W(coordinates));
        denominator = sum(Weights(coordinates));
        
        if denominator ~= 0.0
            u(icp) = u_numerator/denominator;
            v(icp) = v_numerator/denominator;
            w(icp) = w_numerator/denominator;
        else
            u(icp) = NaN;
            v(icp) = NaN;
            w(icp) = NaN;
            warning('no control points influence this (x,y)')
        end
    end
    clear W U V W
    transformedControlPoints = [u,v,w];
    
end

end




