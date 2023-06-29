function taulist = InverseDynamicsThighZero(thetalist, dthetalist, ...
                                   ddthetalist, g, Ftip, Mlist, Glist, ...
                                   Slist)
n = size(thetalist, 1);
Mi = eye(4);
Ai = sym(zeros(6, n));
AdTi = sym(zeros(6, 6, n + 1));
Vi = sym(zeros(6, n + 1));
Vdi = sym(zeros(6, n + 1));
Vdi(4: 6, 1) = -g;
AdTi(:, :, n + 1) = Adjoint(TransInv(Mlist(:, :, n + 1)));
Fi = Ftip;
taulist = sym(zeros(n, 1));
for i=1: n    
    Mi = Mi * Mlist(:, :, i);
    Ai(:, i) = Adjoint(TransInv(Mi)) * Slist(:, i);
    if i == 2
        AdTi(:, :, i) = Adjoint(MatrixExp6(VecTose3(Ai(:, i) ...
                        * -0.0)) * TransInv(Mlist(:, :, i)));
    else
        AdTi(:, :, i) = Adjoint(MatrixExp6Symbolic(VecTose3(Ai(:, i) ...
                        * -thetalist(i)), thetalist(i)) * TransInv(Mlist(:, :, i)));
    end    
    Vi(:, i + 1) = AdTi(:, :, i) * Vi(:, i) + Ai(:, i) * dthetalist(i);
    Vdi(:, i + 1) = AdTi(:, :, i) * Vdi(:, i) ...
                    + Ai(:, i) * ddthetalist(i) ...
                    + ad(Vi(:, i + 1)) * Ai(:, i) * dthetalist(i);    
end
for i = n: -1: 1
    Fi = AdTi(:, :, i + 1)' * Fi + Glist(:, :, i) * Vdi(:, i + 1) ...
         - ad(Vi(:, i + 1))' * (Glist(:, :, i) * Vi(:, i + 1));
    taulist(i) = Fi' * Ai(:, i);
end
end