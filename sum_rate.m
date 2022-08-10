function result = sum_rate(p, xTy, xTz, lambda, alpha)
    p = reshape(p, [3,3,3]);
    prob = zeros(3, 3, 3, 3, 3);
    for u = 1:3
        for v = 1:3
            for x = 1:3
                for y = 1:3
                    for z = 1:3
                        prob(u,v,x,y,z) = p(u,v,x) * xTy(x,y) * xTz(x,z);
                    end
                end
            end
        end
    end
    px = sum(prob,[1,2,4,5]);
    px = squeeze(px);
    %disp('X');
    %disp(px);
    py = sum(prob,[1,2,3,5]);
    py = squeeze(py);
    %disp('Y');
    %disp(py);
    pz = sum(prob,[1,2,3,4]);
    pz = squeeze(pz);
    %disp('Z');
    %disp(pz);
    pu = sum(prob,[2,3,4,5]);
    pu = squeeze(pu);
    %disp('U');
    %disp(pu);
    pxy = sum(prob,[1,2,5]);
    pxy = squeeze(pxy);
    %disp('XY');
    %disp(pxy);
    pxz = sum(prob,[1,2,4]);
    pxz = squeeze(pxz);
    %disp('XZ');
    %disp(pxz);
    puy = sum(prob,[2,3,5]);
    puy = squeeze(puy);
    %disp('UY');
    %disp(puy);
    pvz = sum(prob,[1,3,4]);
    pvz = squeeze(pvz);
    %disp('VZ');
    %disp(pvz);
    puv = sum(prob,[3,4,5]);
    puv = squeeze(puv);
    %disp('UV');
    %disp(puv);
    result = -lambda*h(px)+alpha*h(py)+(1-alpha)*h(pz)+(lambda-1)*h(pu)+(lambda-alpha)*h(pxy)+alpha*h(pxz)-lambda*h(puy)-h(pvz)+h(puv);
end