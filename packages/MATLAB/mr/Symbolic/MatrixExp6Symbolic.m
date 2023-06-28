function T = MatrixExp6Symbolic(se3mat, theta)

omgmat = se3mat(1: 3, 1: 3) / abs(theta); 
T = [MatrixExp3Symbolic(se3mat(1: 3, 1: 3), theta), ...
     (eye(3) * abs(theta) + (1 - cos(abs(theta))) * omgmat ...
      + (abs(theta) - sin(abs(theta))) * omgmat * omgmat) ...
        * se3mat(1: 3, 4) / abs(theta);
     0, 0, 0, 1];
end