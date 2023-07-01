function T = MatrixExp6Symbolic(se3mat, theta)

omgmat = se3mat(1: 3, 1: 3); 
T = [MatrixExp3Symbolic(se3mat(1: 3, 1: 3), theta), ...
     (eye(3) * theta + (1 - cos(theta)) * omgmat ...
      + (theta - sin(theta)) * omgmat * omgmat) ...
        * se3mat(1: 3, 4);
     0, 0, 0, 1];

end