function  R = MatrixExp3Symbolic(so3mat, theta)

omgmat = so3mat / abs(theta);
R = eye(3) + sin(abs(theta)) * omgmat + (1 - cos(abs(theta))) * omgmat * omgmat;
end