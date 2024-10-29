function [u,v] = GetUV_LWLinear_2D(x,y,tdata,W)
%  Given a location in output space (x,y), find the mapped to coordinates
%  (u,v)

T = tdata.tdata.LWLinearTData;
xCP = tdata.tdata.ControlPoints(:,1);
yCP = tdata.tdata.ControlPoints(:,2);
nControlPoints = size(W,1);

u = zeros(nControlPoints,1);
v = zeros(nControlPoints,1);

for icp = 1:nControlPoints
    
    if W(icp) ~= 0.0
        u(icp) = T(1,1,icp)*(x-xCP(icp)) + T(2,1,icp)*(y-yCP(icp)) + T(3,1,icp);
        v(icp) = T(1,2,icp)*(x-xCP(icp)) + T(2,2,icp)*(y-yCP(icp)) + T(3,2,icp);
        
    else
        u(icp) = NaN;
        v(icp) = NaN;
        
    end
    
end
end

