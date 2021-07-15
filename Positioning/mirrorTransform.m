function XyzCoor_Mirror=mirrorTransform(XyzCoordinate,XyzCoordinateCenter)
XyzCoor_Mirror = ones(size(XyzCoordinate));
N = XyzCoordinateCenter*(1./norm(XyzCoordinateCenter));

Trans=[ 2*N(1)*N(1)-1   2*N(1)*N(2)     2*N(1)*N(3); ...
        2*N(1)*N(2)     2*N(2)*N(2)-1   2*N(2)*N(3); ...
        2*N(1)*N(3)     2*N(2)*N(3)     2*N(3)*N(3)-1];
XyzCoor_Mirror = XyzCoordinate*Trans;

% check difference
Check1 = norm(XyzCoordinate-XyzCoordinateCenter);
Check2 = norm(XyzCoor_Mirror-XyzCoordinateCenter);
diff_check=Check2-Check1;
% disp(['映射卫星位置伪距误差： ',num2str(diff_check)]);
end