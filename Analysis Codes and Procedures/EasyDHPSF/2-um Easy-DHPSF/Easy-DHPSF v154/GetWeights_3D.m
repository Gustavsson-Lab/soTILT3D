function [W, coordinates] = GetWeights_3D(x,y,z,tdata,weight_type)
%  This is the matlab code for the C function GetW in inv_lwm.c
%  Given a location in output space (x,y,z), find the weights of each 
%  contributing control point 

xControlPoints = tdata.tdata.ControlPoints(:,1);
yControlPoints = tdata.tdata.ControlPoints(:,2);
zControlPoints = tdata.tdata.ControlPoints(:,3);

dx = x - xControlPoints;
dy = y - yControlPoints;
dz = z - zControlPoints;
dist_to_cp = sqrt( dx.*dx + dy.*dy + dz.*dz);

Ri = dist_to_cp ./ tdata.tdata.RadiusToKthNeighbor;			% This is the 12th closest point by Matlab default

nControlPoints = size(Ri,1);
W = zeros(nControlPoints,1);

% compute the weighting function
switch weight_type
    case 'Maude'
        Ri(Ri>=1) = 0.0;
        [coordinates, ~] = find(Ri);
        Ri2 = Ri(coordinates).*Ri(coordinates);
        Ri3 = Ri(coordinates).*Ri2;
        W(coordinates) = 1.0 - 3.0*Ri2 + 2.0*Ri3;           % weight of ControlPoint i
    case 'Gaussian'                                         % see Goshtasby, 2012, Eqn 9.91
        sigma =  tdata.tdata.RadiusToKthNeighbor;               
        s =  tdata.tdata.smoothnessParam;                   % smoothness parameter 
        W = exp(-( (dx.^2 + dy.^2+ dz.^2) ./ (2*(s*sigma).^2) ));
        coordinates = 1:nControlPoints;
    otherwise
        error(message('GetWeights:wrong weight type'))
end

end

