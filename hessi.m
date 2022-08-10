% Hessian

function [Da, Db, D] = hessi(p1, p2, f1, f2, xTy1, xTy2, xTz1, xTz2, lambda, alpha)
    %p = reshape(p, [9, 9, 9]);
    p1 = reshape(p1, [2, 2]);
    p2 = reshape(p2, [2, 2]);
    %disp('p1')
    %disp(p1)
    %disp('p2')
    %disp(p2)
    p = kron(p1, p2);
    %disp('p')
    %disp(p)
    %disp('xTy1')
    %disp(xTy1)
    %disp('xTy2')
    %disp(xTy2)
    %disp('xTz1')
    %disp(xTz1)
    %disp('xTz2')
    %disp(xTz2)
    prob1 = zeros(2, 2, 3, 3, 3);
    prob2 = zeros(2, 2, 3, 3, 3);
    for u = 1:2
        for v = 1:2
            x1 = f1(u, v);
            x2 = f2(u, v);
            for x = 1:3
                for y = 1:3
                    for z = 1:3
                        if x == x1 
                            prob1(u, v, x, y, z) = p1(u, v) * xTy1(x, y) * xTz1(x, z);
                        end
                        if x == x2
                            prob2(u, v, x, y, z) = p2(u, v) * xTy2(x, y) * xTz2(x, z);
                        end 
                    end
                end
            end
        end
    end
    sizeA = 4;
    Auv1 = zeros(sizeA, sizeA);
    Auv2 = zeros(sizeA, sizeA);
    Auy1 = zeros(sizeA, sizeA);
    Auy2 = zeros(sizeA, sizeA);
    Avz1 = zeros(sizeA, sizeA);
    Avz2 = zeros(sizeA, sizeA);
    Ay1 = zeros(sizeA, sizeA);
    Ay2 = zeros(sizeA, sizeA);
    Az1 = zeros(sizeA, sizeA);
    Az2 = zeros(sizeA, sizeA);
    Au1 = zeros(sizeA, sizeA);
    Au2 = zeros(sizeA, sizeA);
    p1 = reshape(p1, [1,4]);
    p2 = reshape(p2, [1,4]);
    for i = 1:sizeA
        Auv1(i, i) = p1(i);
        Auv2(i, i) = p2(i);
    end
    %disp('Auv1')
    %disp(Auv1)
    %disp('Auv2')
    %disp(Auv2)
    puvy1 = sum(prob1, [3, 5]);
    puvy1 = squeeze(puvy1);
    puvy2 = sum(prob2, [3, 5]);
    puvy2 = squeeze(puvy2);
    puy1 = sum(prob1, [2, 3, 5]);
    puy1 = squeeze(puy1);
    puy2 = sum(prob2, [2, 3, 5]);
    puy2 = squeeze(puy2);
    for i = 1:sizeA
        for j = 1:sizeA
            [u1, v1] = ind2sub([2, 2], i);
            [u2, v2] = ind2sub([2, 2], j);
            for y = 1:3
                if u1 == u2 && puy1(u1, y) > 1e-5
                    Auy1(i, j) = Auy1(i, j) + puvy1(u1, v1, y) * puvy1(u2, v2, y) / puy1(u1, y);
                end
                if u1 == u2 && puy2(u2, y) > 1e-5
                    Auy2(i, j) = Auy2(i, j) + puvy2(u1, v1, y) * puvy2(u2, v2, y) / puy2(u1, y);
                end
            end
        end
    end
    %disp('Auy1')
    %disp(Auy1)
    %disp('Auy2')
    %disp(Auy2)
    puvz1 = sum(prob1, [3, 4]);
    puvz1 = squeeze(puvz1);
    puvz2 = sum(prob2, [3, 4]);
    puvz2 = squeeze(puvz2);
    pvz1 = sum(prob1, [1, 3, 4]);
    pvz1 = squeeze(pvz1);
    pvz2 = sum(prob2, [1, 3, 4]);
    pvz2 = squeeze(pvz2);
    for i = 1:sizeA
        for j = 1:sizeA
            [u1, v1] = ind2sub([2, 2], i);
            [u2, v2] = ind2sub([2, 2], j);
            for z = 1:3
                 if v1 == v2 && pvz1(v1, z) > 1e-5
                     Avz1(i, j) = Avz1(i, j) + puvz1(u1, v1, z) * puvz1(u2, v2, z) / pvz1(v1, z);
                 end
                 if v1 == v2 && pvz2(v1, z) > 1e-5
                     Avz2(i, j) = Avz2(i, j) + puvz2(u1, v1, z) * puvz2(u2, v2, z) / pvz2(v1, z);
                 end
            end
        end
    end
    %disp('Avz1')
    %disp(Avz1)
    %disp('Avz2')
    %disp(Avz2)
    py1 = sum(prob1, [1, 2, 3, 5]);
    py1 = squeeze(py1);
    py2 = sum(prob2, [1, 2, 3, 5]);
    py2 = squeeze(py2);
    for i = 1:sizeA
        for j = 1:sizeA
            [u1, v1] = ind2sub([2, 2], i);
            [u2, v2] = ind2sub([2, 2], j);
            for y = 1:3
                if py1(y) > 1e-5
                    Ay1(i, j) = Ay1(i, j) + puvy1(u1, v1, y) * puvy1(u2, v2, y) / py1(y);
                end
                if py2(y) > 1e-5
                    Ay2(i, j) = Ay2(i, j) + puvy2(u1, v1, y) * puvy2(u2, v2, y) / py2(y);
                end
            end
        end
    end
    %disp('Ay1')
    %disp(Ay1)
    %disp('Ay2')
    %disp(Ay2)
    pz1 = sum(prob1, [1, 2, 3, 4]);
    pz1 = squeeze(pz1);
    pz2 = sum(prob2, [1, 2, 3, 4]);
    pz2 = squeeze(pz2);
    for i = 1:sizeA
        for j = 1:sizeA
            [u1, v1] = ind2sub([2, 2], i);
            [u2, v2] = ind2sub([2, 2], j);
            for z = 1:3
               if pz1(z) > 1e-5
                   Az1(i, j) = Az1(i, j) + puvz1(u1, v1, z) * puvz1(u2, v2, z) / pz1(z);
               end
               if pz2(z) > 1e-5
                   Az2(i, j) = Az2(i, j) + puvz2(u1, v1, z) * puvz2(u2, v2, z) / pz2(z);
               end
            end
        end
    end
    %disp('Az1')
    %disp(Az1)
    %disp('Az2')
    %disp(Az2)
    pu1 = sum(prob1, [2, 3, 4, 5]);
    pu1 = squeeze(pu1);
    pu2 = sum(prob2, [2, 3, 4, 5]);
    pu2 = squeeze(pu2);
    puv1 = sum(prob1, [3, 4, 5]);
    puv1 = squeeze(puv1);
    puv2 = sum(prob2, [3, 4, 5]);
    puv2 = squeeze(puv2);
    for i = 1:sizeA
        for j = 1:sizeA
            [u1, v1] = ind2sub([2, 2], i);
            [u2, v2] = ind2sub([2, 2], j);
            if u1 == u2 && pu1(u) > 1e-5
                Au1(i, j) = Au1(i, j) + puv1(u1, v1) * puv1(u2, v2) / pu1(u1);
            end
            if u1 == u2 && pu2(u) > 1e-5
                Au2(i, j) = Au2(i, j) + puv2(u1, v1) * puv2(u2, v2) / pu2(u1);
            end
        end
    end
    %disp('Au1')
    %disp(Au1)
    %disp('Au2')
    %disp(Au2)
    Auv = zeros(sizeA * sizeA, sizeA * sizeA);
    Au = zeros(sizeA * sizeA, sizeA * sizeA);
    Ay = zeros(sizeA * sizeA, sizeA * sizeA);
    Az = zeros(sizeA * sizeA, sizeA * sizeA);
    Auy = zeros(sizeA * sizeA, sizeA * sizeA);
    Avz = zeros(sizeA * sizeA, sizeA * sizeA);
    Auv = kron(Auv1, Auv2);
    Au = kron(Au1, Au2);
    Ay = kron(Ay1, Ay2);
    Az = kron(Az1, Az2);
    Auy = kron(Auy1, Auy2);
    Avz = kron(Avz1, Avz2);
    A = zeros(sizeA * sizeA, sizeA * sizeA);
    B = zeros(sizeA * sizeA, sizeA * sizeA);
    hess = zeros(sizeA * sizeA, sizeA * sizeA);
    A = Auv + (lambda - 1) * Au + alpha * Ay + (1 - alpha) * Az;
    B = Auy + Avz;
    hess = A - B;
    [Va, Da] = eig(A);
    [Vb, Db] = eig(B);
    [V, D] = eig(hess);
    %disp('diagonal of A')
    %disp(Da)
    %disp('diagonal of B')
    %disp(Db)
    %disp('diagonal of hessian')
    %disp(D)
    %puv = sum(p, 3);
    %puv = squeeze(puv);
    %eps = 1e-7;
    %index = find(sum(puv, 1) < eps);
    %puv(:, index) = [];
    %p(:, index, :) = [];
    %index = find(sum(puv, 2) < eps);
    %puv(index, :) = [];
    %p(index, :, :) = [];
    %size1 = size(p, 1);
    %size2 = size(p, 2);
    %size3 = size(p,3);
    %disp('size u')
    %disp(size1)
    %disp('size v')
    %disp(size2)
    %disp('size x')
    %disp(size3)
    %for u = 1:size1
    %    for v = 1:size2
    %        for x = 1:size3
    %            for y = 1:3
    %                for z = 1:3
    %                    prob(u,v,x,y,z) = p(u,v,x) * xTy(x,y) * xTz(x,z);
    %                    disp('   [ u,    v,    x,    y,    z]')
    %                    disp([u, v, x, y, z])
    %                end
    %            end
    %        end
    %    end
    %end
    %px = sum(prob,[1,2,4,5]);
    %px = squeeze(px);
    %disp('X');
    %disp(px);
    %py = sum(prob,[1,2,3,5]);
    %py = squeeze(py);
    %disp('Y');
    %disp(py);
    %pz = sum(prob,[1,2,3,4]);
    %pz = squeeze(pz);
    %disp('Z');
    %disp(pz);
    %pu = sum(prob,[2,3,4,5]);
    %pu = squeeze(pu);
    %disp('U');
    %disp(pu);
    %pxy = sum(prob,[1,2,5]);
    %pxy = squeeze(pxy);
    %disp('XY');
    %disp(pxy);
    %pxz = sum(prob,[1,2,4]);
    %pxz = squeeze(pxz);
    %disp('XZ');
    %disp(pxz);
    %puy = sum(prob,[2,3,5]);
    %puy = squeeze(puy);
    %disp('UY');
    %disp(puy);
    %pvz = sum(prob,[1,3,4]);
    %pvz = squeeze(pvz);
    %disp('VZ');
    %disp(pvz);
    %puv = sum(prob,[3,4,5]);
    %puv = squeeze(puv);
    %disp('UV');
    %disp(puv);
    %puvxy = sum(prob, 5);
    %puvxy = squeeze(puvxy);
    %disp('UVXY')
    %disp(puvxy)
    %puvxz = sum(prob, 4);
    %puvxz = squeeze(puvxz);
    %disp('UVXZ')
    %disp(puvxz)
    %puvx = sum(prob, [4,5]);
    %puvx = squeeze(puvx);
    %disp('UVX');
    %disp(puvx);
    %Auv = zeros(27,27);
    %disp('size of u')
    %disp(size1)
    %disp('size of v')
    %disp(size2)
    %sizeA = size1 * size2;
    %Auv = zeros(sizeA, sizeA);
    %for i = 1:sizeA
    %    for j = 1:sizeA
    %        [u1, v1, x1] = ind2sub([3, 3, 3], i);
    %        [u2, v2, x2] = ind2sub([3, 3, 3], j);
    %        for y = 1:3
    %            if u1 == u2 && v1 == v2 && puv(u1, v1) > 1e-5
    %                Auv(i, j) = Auv(i, j) + puvx(u1, v1, x1) * puvx(u2, v2, x2) / puv(u1, v1);
    %            end
    %        end
    %    end
    %end
    %disp('Auv')
    %disp(Auv)
    %Auy = zeros(sizeA,sizeA);
    %for i = 1:sizeA
    %    for j = 1:sizeA
    %        [u1, v1, x1] = ind2sub([3, 3, 3], i);
    %        [u2, v2, x2] = ind2sub([3, 3, 3], j);
    %        for y = 1:3
    %            if  u1 == u2 && puy(u1, y) > 1e-5
    %                Auy(i, j) = Auy(i, j) + puvxy(u1, v1, x1, y) * puvxy(u2, v2, x2, y) / puy(u1, y);
    %            end
    %        end
    %    end
    %end
    %disp('Auy')
    %disp(Auy)
    %Avz = zeros(sizeA,sizeA);
    %for i = 1:sizeA
    %    for j = 1:sizeA
    %        [u1, v1, x1] = ind2sub([3, 3, 3], i);
    %        [u2, v2, x2] = ind2sub([3, 3, 3], j);
    %        for z = 1:3
    %            if v1 == v2 && pvz(v1, z) > 1e-5
    %                Avz(i, j) = Avz(i, j) + puvxz(u1, v1, x1, z) * puvxz(u2, v2, x2,z) / pvz(v1, z);
    %            end
    %        end
    %    end
    %end
    %disp('Avz')
    %disp(Avz)    
    %Ay = zeros(sizeA, sizeA);
    %for i = 1:sizeA
    %    for j = 1:sizeA 
    %        [u1, v1, x1] = ind2sub([3, 3, 3], i);
    %        [u2, v2, x2] = ind2sub([3, 3, 3], j);
    %        for y = 1:3
    %            if py(y) > 1e-5
    %                Ay(i, j) = Ay(i, j) + puvxy(u1, v1, x1, y) * puvxy(u2, v2, x2, y) / py(y);
    %            end
    %        end
    %    end
    %end
    %disp('Ay')
    %disp(Ay)
    %Az = zeros(sizeA, sizeA);
    %for i = 1:sizeA
    %    for j = 1:sizeA
    %        [u1, v1, x1] = ind2sub([3, 3, 3], i);
    %        [u2, v2, x2] = ind2sub([3, 3, 3], j);
    %        for z = 1:3
    %            if pz(z) > 1e-5
    %                Az(i, j) = Az(i, j) + puvxz(u1, v1, x1, z) * puvxz(u2, v2, x2, z) / pz(z);
    %            end
    %        end
    %    end
    %end
    %disp('Az')
    %disp(Az)
    %Au = zeros(sizeA, sizeA);
    %for i = 1:sizeA
    %    for j = 1:sizeA
    %        [u1, v1, x1] = ind2sub([3, 3, 3], i);
    %        [u2, v2, x2] = ind2sub([3, 3, 3], j);
    %        if u1 == u2 && pu(u) > 1e-5
    %            Au(i, j) = Au(i, j) + puvx(u1, v1, x1) * puvx(u2, v2, x2) / pu(u);
    %        end
    %    end
    %end
    %disp('Au')
    %disp(Au)
    
    %hess = lambda * Auy + Avz - alpha * Ay - (1-alpha) * Az - (lambda - 1) * Au - Auv;
    %disp('hess')
    %disp(hess)
    %[V, D] = eig(hess);
    %disp('Diagonal')
    %disp(D)
end