function [u,v,w] = GetUV_LWQuadratic_3D(x,y,z,tdata,W)
%  Given a location in output space (x,y,z), find the mapped to coordinates
%  (u,v,w)

T = tdata.tdata.LWQuadraticTData;
xCP = tdata.tdata.ControlPoints(:,1);
yCP = tdata.tdata.ControlPoints(:,2);
zCP = tdata.tdata.ControlPoints(:,3);
nControlPoints = size(W,1);

u = zeros(nControlPoints,1);
v = zeros(nControlPoints,1);
w = zeros(nControlPoints,1);

for icp = 1:nControlPoints
    
    if W(icp) ~= 0.0
        
        u(icp) = T(1,1,icp)*(x-xCP(icp)) + T(2,1,icp)*(y-yCP(icp)) + T(3,1,icp)*(z-zCP(icp)) + T(4,1,icp)*(x-xCP(icp))^2 + T(5,1,icp)*(y-yCP(icp))^2 + T(6,1,icp)*(z-zCP(icp))^2 + T(7,1,icp)*(x-xCP(icp))*(y-yCP(icp)) + T(8,1,icp)*(x-xCP(icp))*(z-zCP(icp)) + T(9,1,icp)*(y-yCP(icp))*(z-zCP(icp)) + T(10,1,icp);
        v(icp) = T(1,2,icp)*(x-xCP(icp)) + T(2,2,icp)*(y-yCP(icp)) + T(3,2,icp)*(z-zCP(icp)) + T(4,2,icp)*(x-xCP(icp))^2 + T(5,2,icp)*(y-yCP(icp))^2 + T(6,2,icp)*(z-zCP(icp))^2 + T(7,2,icp)*(x-xCP(icp))*(y-yCP(icp)) + T(8,2,icp)*(x-xCP(icp))*(z-zCP(icp)) + T(9,2,icp)*(y-yCP(icp))*(z-zCP(icp)) + T(10,2,icp);
        w(icp) = T(1,3,icp)*(x-xCP(icp)) + T(2,3,icp)*(y-yCP(icp)) + T(3,3,icp)*(z-zCP(icp)) + T(4,3,icp)*(x-xCP(icp))^2 + T(5,3,icp)*(y-yCP(icp))^2 + T(6,3,icp)*(z-zCP(icp))^2 + T(7,3,icp)*(x-xCP(icp))*(y-yCP(icp)) + T(8,3,icp)*(x-xCP(icp))*(z-zCP(icp)) + T(9,3,icp)*(y-yCP(icp))*(z-zCP(icp)) + T(10,3,icp);
        
    else
        
        u(icp) = NaN;
        v(icp) = NaN;
        w(icp) = NaN;
        
    end
    
end
end

