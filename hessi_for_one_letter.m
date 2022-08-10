function hessi_for_one_letter(p, fuv, xTy, xTz, lambda, alpha)
    p = reshape(p, [2, 2])
    prob = zeros(2, 2, 3, 3, 3);
    for u = 1:2
        for v = 1:2
            x0 = fuv(u, v);
            for x = 1:3
                for y = 1:3
                    for z = 1:3
                        if x == x0 
                            prob(u, v, x, y, z) = p(u, v) * xTy(x, y) * xTz(x, z);
                        end
                    end
                end
            end
        end
    end
    Auv = zeros(4, 4);
    Auy = zeros(4, 4);
    Avz = zeros(4 ,4);
    Au = zeros(4, 4);
    Ay = zeros(4, 4);
    Az = zeros(4, 4);
    p0 = reshape(p, [1, 4]);
    sizeA = 4;
    for i = 1:sizeA
        Auv(i, i) = p0(i);
    end
    puvy = sum(prob, [3, 5]);
    puvy = squeeze(puvy);
    puy = sum(prob, [2, 3, 5]);
    puy = squeeze(puy);
    for i = 1:sizeA
        for j = 1:sizeA
            [u1, v1] = ind2sub([2, 2], i);
            [u2, v2] = ind2sub([2, 2], j);
            for y = 1:3
                if u1 == u2 && puy(u1, y) > 1e-5
                    Auy(i, j) = Auy(i, j) + puvy(u1, v1, y) * puvy(u2, v2, y) / puy(u1, y);
                end
            end
        end
    end
    puvz = sum(prob, [3, 4]);
    puvz = squeeze(puvz);
    pvz = sum(prob, [1, 3, 4]);
    pvz = squeeze(pvz);
    for i = 1:sizeA
        for j = 1:sizeA
            [u1, v1] = ind2sub([2, 2], i);
            [u2, v2] = ind2sub([2, 2], j);
            for z = 1:3
                 if v1 == v2 && pvz(v1, z) > 1e-5
                     Avz(i, j) = Avz(i, j) + puvz(u1, v1, z) * puvz(u2, v2, z) / pvz(v1, z);
                 end
            end
        end
    end
    py = sum(prob, [1, 2, 3, 5]);
    py = squeeze(py);
    for i = 1:sizeA
        for j = 1:sizeA
            [u1, v1] = ind2sub([2, 2], i);
            [u2, v2] = ind2sub([2, 2], j);
            for y = 1:3
                if py(y) > 1e-5
                    Ay(i, j) = Ay(i, j) + puvy(u1, v1, y) * puvy(u2, v2, y) / py(y);
                end
            end
        end
    end
    pz = sum(prob, [1, 2, 3, 4]);
    pz = squeeze(pz);
    for i = 1:sizeA
        for j = 1:sizeA
            [u1, v1] = ind2sub([2, 2], i);
            [u2, v2] = ind2sub([2, 2], j);
            for z = 1:3
               if pz(z) > 1e-5
                   Az(i, j) = Az(i, j) + puvz(u1, v1, z) * puvz(u2, v2, z) / pz(z);
               end
            end
        end
    end
    pu = sum(prob, [2, 3, 4, 5]);
    pu = squeeze(pu);
    puv = sum(prob, [3, 4, 5]);
    puv = squeeze(puv);
    for i = 1:sizeA
        for j = 1:sizeA
            [u1, v1] = ind2sub([2, 2], i);
            [u2, v2] = ind2sub([2, 2], j);
            if u1 == u2 && pu(u) > 1e-5
                Au(i, j) = Au(i, j) + puv(u1, v1) * puv(u2, v2) / pu(u1);
            end
        end
    end
    A = zeros(sizeA * sizeA, sizeA * sizeA);
    B = zeros(sizeA * sizeA, sizeA * sizeA);
    hess = zeros(sizeA * sizeA, sizeA * sizeA);
    A = Auv + (lambda - 1) * Au + alpha * Ay + (1 - alpha) * Az;
    B = Auy + Avz;
    hess = A - B;
    [Va, Da] = eig(A);
    [Vb, Db] = eig(B);
    [V, D] = eig(hess);
    disp('diagonal of A')
    disp(Da)
    disp('diagonal of B')
    disp(Db)
    disp('diagonal of hessian')
    disp(D)
end