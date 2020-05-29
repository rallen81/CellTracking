function Rout = rotate(X,Y,rotate)
% rotate by 'rotate' + is anticlockwise)

R=(X.^2+Y.^2).^(0.5);

theta=atan2(Y,X);

theta=theta+rotate;

Xrot=R.*cos(theta);
Yrot=R.*sin(theta);
Rout=[Xrot Yrot];
end

