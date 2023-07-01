function T = FKinSpaceSymbolic(M, Slist, thetalist)

T = M;
for i = size(thetalist): -1: 1
    T = MatrixExp6Symbolic(VecTose3(Slist(:, i)), thetalist(i)) * T;
end
end