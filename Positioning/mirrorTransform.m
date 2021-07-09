function XyzCoor_Mirror=mirrorTransform(XyzCoordinate1,XyzCoordinate2)
XyzCoor_Mirror=ones(size(XyzCoordinate1));
N=XyzCoordinate2*(1./norm(XyzCoordinate2));
Trans=[2*N(1)*N(1)-1 2*N(1)*N(2) 2*N(1)*N(3); 2*N(1)*N(2) 2*N(2)*N(2)-1 2*N(2)*N(3); 2*N(1)*N(3) 2*N(2)*N(3) 2*N(3)*N(3)-1]
XyzCoor_Mirror=XyzCoordinate1*Trans;
Check1=norm(XyzCoordinate1-XyzCoordinate2);
Check2=norm(XyzCoor_Mirror-XyzCoordinate2);
diff_check=Check2-Check1;
disp(['映射卫星位置伪距误差： ',num2str(diff_check)]);
end