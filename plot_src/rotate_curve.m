function [yp, zp] = rotate_curve(y,z,theta)

% fprintf("y = %3.3f and z = %3.3f\n",y, z);

yp = cos(theta)*y - sin(theta)*z;
zp = sin(theta)*y + cos(theta)*z;
% fprintf("y = %3.3f and z = %3.3f\n",yp, zp);
