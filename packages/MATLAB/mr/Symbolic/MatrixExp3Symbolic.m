function  R = MatrixExp3Symbolic(so3mat, theta)

omgmat = so3mat;
R = eye(3) + sin(theta) * omgmat + (1 - cos(theta)) * omgmat * omgmat;
end