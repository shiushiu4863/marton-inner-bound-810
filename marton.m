clear;
close all;

m = 3;
for amount = 1:10
    lambda = 2 * rand + 1;
    flag = 0;
    alpha = rand;
    filename = sprintf('%.2foutput%.2f.txt', lambda, alpha);
    output_file = fopen(filename, 'a+');
    sample = 1;
    for trial=1:100
        if flag == 1
            break
        end
        
        xTy = random_channel(m, m);
        xTz = random_channel(m, m);

        xTy2 = kron(xTy, xTy);
        xTz2 = kron(xTz, xTz);
    
        N = 3 * 3 * 3;
        options = optimoptions('fmincon','algorithm','sqp','TolFun',1e-12,'TolCon',1e-12,'TolX',1e-12,'MaxIter',1e12,'Display','off'); 

        eps = 0;
        ms = MultiStart('FunctionTolerance',1e-8,'UseParallel',true);
        x0 = ones(1, N) / N;
        problem = createOptimProblem('fmincon','x0',x0,...
            'objective', @(x) -sum_rate(x, xTy, xTz, lambda, alpha) ,...
            'Aeq', ones(1, N), 'beq', 1, 'lb', zeros(1, N), 'ub', ones(1, N), 'options',options);
        numParallel = 100;
        [x, fval, eflag, output, manymins] = run(ms, problem, numParallel);
        for i = 1:length(manymins)
            if flag == 1
                break
            end
            solution = manymins(i);

            if -solution.Fval > 0
                disp(solution);
                disp("Channel W(Y | X):")
                disp(xTy);
            
                disp("Channel W(Z | X):")
                disp(xTz);
                X = solution.X;            
                %disp("Proposed local maxima:")
                %disp(X);
                %disp("Proposed local maxima:(array)")
                %disp(reshape(X, [2,2,3]));
                index = find(sum(X, 1) < 1e-6);
                X(:, index, :) = [];
                index = find(sum(X, 2) < 1e-6);
                X(index, :, :) = [];
                %disp('X')
                %disp(X(:)')
                %disp('check')
                %disp(reshape(X, [1,27]));
                %disp(grad(solution.X, xTy, xTz, lambda, alpha))
                [G, result2] = grad(solution.X, xTy, xTz, lambda, alpha);
                G = G(:)';
                %disp('Gradient')
                %disp(G)
                %p_opt = reshape(solution.X, [3, 3, 3]);
                %pTp = kron(p_opt, p_opt)
                pTp = kron(solution.X, solution.X);
                %disp('p*')
                %disp(solution.X)
                solution = reshape(solution.X, [3, 3, 3]);
                solution = sum(solution, 3);
                %disp('size(pTp)')
                %disp(size(pTp))
                %disp('pTp')
                %disp(pTp)
                %hessian = hessi(pTp, xTy, xTz, lambda, alpha);
                %disp('hess')
                %disp(hessian)
                %hessian = hess(solution.X, xTy2, xTz2, lambda);
                %hessian_trunc = hess_trunc(X, hessian);
                %disp('Hessian');
                %disp(hessian);
                %disp('Hessian Truncated');
                %disp(hessian_trunc);
                %disp('Eig:')
                %disp(eig(hessian_trunc));
                %if max(eig(hessian_trunc), [], 'all') <= 0
                %    flag = 1;
                %    break;
                %end
            
                epsilon = 1e-7;
                index = find(sum(solution, 1) < epsilon);
                solution(:, index) = [];
                index = find(sum(solution, 2) < epsilon);
                solution(index, :) = [];
                %disp('solution')
                %disp(solution2)
                size1 = size(solution, 1);
                size2 = size(solution, 2);
                size3 = size(solution, 3);
                pxuv = zeros(size3, size1, size2);
                %for x = 1:size3
                %    for u = 1:size1
                %        for v = 1:size2
                %            pxuv(x, u, v) = solution(u, v, x) / solution(u, v);
                %         end
                %   end
                %end
                %disp('px|uv')
                %disp(pxuv)
                disp('result')
                disp(result2)
                if size(solution) == [2, 2]
                    fuv = zeros(2, 2);
                    for x = 1:3 
                        for u = 1:2
                            for v = 1:2
                                if result2(x, u, v) > 0.9 
                                    fuv(u, v) = x;
                                end
                            end
                        end
                    end
                    fprintf(output_file, 'Trial %d\n', sample);
                    sample = sample + 1;
                    fprintf(output_file, '%12s\n', 'xTy');
                    fprintf(output_file, '%.4f %.4f %.4f\n%.4f %.4f %.4f\n%.4f %.4f %.4f\n', xTy);
                    fprintf(output_file, '%12s\n', 'xTz');
                    fprintf(output_file, '%.4f %.4f %.4f\n%.4f %.4f %.4f\n%.4f %.4f %.4f\n', xTz);
                    fprintf(output_file, '%6s\n', 'p*');
                    fprintf(output_file, '%.4f %.4f\n%.4f %.4f\n', solution);
                    fprintf(output_file, '%6s\n', 'fuv');
                    fprintf(output_file, '%.4f %.4f\n%.4f %.4f\n', fuv);
                    fprintf(output_file, '%2s\n', 'lambda');
                    fprintf(output_file, '%.2f\n', lambda);
                    fprintf(output_file, '%2s\n', 'alpha');
                    fprintf(output_file, '%.2f\n\n', alpha);
                    disp('  file  trial')
                    disp([amount, trial])
                    break;
                end 
            end
        end
    end
end

function [result1, result2] = grad(p, xTy, xTz, lambda, alpha)
    p = reshape(p, [3,3,3]);
    %disp('p')
    %disp(p)
    eps = 1e-7;
    prob = zeros(3, 3, 3, 3, 3);
    puv = sum(p, 3);
    %puv = sum(prob, [3, 4, 5]);
    puv = squeeze(puv);
    index = find(sum(puv, 1) < eps);
    puv(:, index) = [];
    p(:, index, :) = [];
    index = find(sum(puv, 2) < eps);
    puv(index, :) = [];
    p(index, : ,:) = [];
    size1 = size(p, 1);
    size2 = size(p, 2);
    size3 = size(p, 3);
    %disp('puv')
    %disp(puv)
    for u = 1:size1
        for v = 1:size2
            for x = 1:size3
                for y = 1:3
                    for z = 1:3
                        prob(u,v,x,y,z) = p(u,v,x) * xTy(x,y) * xTz(x,z);
                    end
                end
            end
        end
    end
    %disp('prob (before)')
    %disp(prob)
    %eps = 1e-4;
    %disp('prob (after)')
    %disp(prob)
    puvx = sum(prob, [4,5]);
    puvx = squeeze(puvx);
    px = sum(prob,[1,2,4,5]);
    px = squeeze(px);
    %disp('px')
    %disp(px)
    py = sum(prob,[1,2,3,5]);
    py = squeeze(py);
    %disp('py')
    %disp(py)
    pz = sum(prob,[1,2,3,4]);
    pz = squeeze(pz);
    %disp('pz')
    %disp(pz)
    pu = sum(prob,[2,3,4,5]);
    pu = squeeze(pu);
    %disp('pu')
    %disp(pu)
    pxy = sum(prob,[1,2,5]);
    pxy = squeeze(pxy);
    %disp('pxy')
    %disp(pxy)
    pxz = sum(prob,[1,2,4]);
    pxz = squeeze(pxz);
    %disp('pxz')
    %disp(pxz)
    puy = sum(prob,[2,3,5]);
    puy = squeeze(puy);
    %disp('puy')
    %disp(puy)
    pvz = sum(prob,[1,3,4]);
    pvz = squeeze(pvz);
    %disp('pvz')
    %disp(pvz)
    puv = sum(prob,[3,4,5]);
    puv = squeeze(puv);
    %disp('puv (before)')
    %disp(puv)
    %eps = 1e-6;
    %index = find(sum(puv,1) < eps);
    %prob(:, index, :, : ,:) = [];
    %puv(:, index) = [];
    %index = find(sum(puv,2) < eps);
    %prob(index, :, :, :, :) = [];
    %puv(index, :) = [];
    %disp('puv (after)')
    %disp(puv)
    result1 = zeros(size(prob, 1:3));
    sizeu = size(prob,1);
    sizev = size(prob,2);
    sizex = size(prob,3);
    sizey = size(prob,4);
    sizez = size(prob,5);
    %disp(size(prob, 1))
    %disp(size(prob, 2))
    %disp(size(prob, 3))
    %disp(size(prob, 4))
    %disp(size(prob, 5))
    %disp('size prob')
    %disp(size(prob))
    for u = 1:sizeu
        for v = 1:sizev
            for x = 1:sizex
                for y = 1:sizey
                    for z = 1:sizez
                        summation = lambda * log2(px(x)) + lambda * log2(puy(u,y)) + log2(pvz(v,z)) - alpha * log2(py(y)) - (1 - alpha) * log2(pz(z)) - (lambda - 1) * log2(pu(u)) - (lambda - alpha) * log2(pxy(x,y)) - alpha * log2(pxz(x,z)) - log2(puv(u,v));
                        %disp('u,v,x,y,z')
                        %disp(u)
                        %disp(v)
                        %disp(x)
                        %disp(y)
                        %disp(z)
                        if ~isnan(summation) && ~isinf(summation)    
                            result1(u,v,x) = result1(u,v,x) + (xTy(x,y) * xTz(x,z) * summation);
                        end
                    end
                end
            end
        end
    end
    %index = find(sum(puv,1) < eps);
    %puv(:, index) = [];
    %index = find(sum(puv,2) < eps);
    %puv(index, :) = [];
    %disp('puv (after)')
    %disp(puv)
    %prob(:, index, :, :, :) = [];
    %prob(index, :, :, :, :) = [];
    pxuv = zeros(size3, size1, size2);
    for x = 1:size3
        for u = 1:size1
            for v = 1:size2
                pxuv(x, u, v) = p(u, v, x) / puv(u, v);
            end
        end
    end
    %disp('px|uv')
    %disp(pxuv)
    result2 = pxuv;
    %disp('detail')
    %disp('px')
    %disp(px)
    %disp('py')
    %disp(py)
    %disp('pz')
    %disp(pz)
    %disp('pu')
    %disp(pu)
    %disp('pxy')
    %disp(pxy)
    %disp('pxz')
    %disp(pxz)
    %disp('puv')
    %disp(puv)
end

function result = is_product(p_X)
    N = sqrt(length(p_X));
    p_X = reshape(p_X, N, N);
    max_error = max(abs(p_X - kron(sum(p_X, 1), sum(p_X, 2))), [], 'all');
    if max_error > 1e-4
        result = 0;
    else
        result = 1;
    end
end


function result = trunc(array, eps)
    for i = 1:size(array, 1)
        for j = 1:size(array, 2)
            if abs(array(i, j)) < eps
                array(i, j) = 0;
            end
        end
    end
    result = array;
end

function result = hess_trunc(p, hess)
    for i = length(p):-1:1
        if p(i) == 0
            hess(i, :) = [];
            hess(:, i) = [];
        end
    end
    result = hess;
end