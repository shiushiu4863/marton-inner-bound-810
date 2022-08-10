function check_hessian()
    input_file = fopen('2.09output0.04.txt', 'r');
    dtb = zeros(100, 26);
    for number = 1:100
        str = fscanf(input_file, '%s', 1);
        %disp(str)
        trial = fscanf(input_file, '%d', 1);
        %disp(trial)
        str = fscanf(input_file, '%s', 1);
        %disp(str)
        xTy = fscanf(input_file, '%f', 9);
        %disp(xTy);
        for i = 1:9
            dtb(number, i) = xTy(i);
        end
        str = fscanf(input_file, '%s', 1);
        %disp(str)
        xTz = fscanf(input_file, '%f', 9);
        %disp(xTz)
        for i = 10:18
            dtb(number, i) = xTz(i - 9);
        end
        str = fscanf(input_file, '%s', 1);
        %disp(str)
        p = fscanf(input_file, '%f', 4);
        %disp(p)
        for i = 19:22
            dtb(number, i) = p(i - 18);
        end
        str = fscanf(input_file, '%s', 1);
        %disp(str)
        fuv = fscanf(input_file, '%f', 4);
        %disp(fuv)
        for i = 23:26
            dtb(number, i) = fuv(i - 22);
        end
        str = fscanf(input_file, '%s', 1);
        %disp(str)
        lambda = fscanf(input_file, '%f', 1);
        %disp(lambda)
        %dtb(number, 27) = lambda;
        str = fscanf(input_file, '%s', 1);
        %disp(str)
        alpha = fscanf(input_file, '%f', 1);
        %disp(alpha)
        %dtb(number, 28) = alpha;
    end
    %disp(dtb)
    %disp('lambda')
    %disp(lambda)
    %disp('alpha')
    %disp(alpha)
    r1 = randi([1, 100]);
    r2 = randi([1, 100]);
    disp('r1')
    disp(r1)
    disp('r2')
    disp(r2)
    r1 = 70;
    r2 = 44;
    p1 = zeros(4, 1);
    p2 = zeros(4, 1);
    fuv1 = zeros(4, 1);
    fuv2 = zeros(4, 1);
    xTy1 = zeros(9, 1);
    xTy2 = zeros(9, 1);
    xTz1 = zeros(9, 1);
    xTz2 = zeros(9, 1);
    for i = 1:9
        xTy1(i) = dtb(r1, i);
        xTy2(i) = dtb(r2, i);
    end
    for i = 10:18
        xTz1(i - 9) = dtb(r1, i);
        xTz2(i - 9) = dtb(r2, i);
    end
    for i = 19:22
        p1(i - 18) = dtb(r1, i);
        p2(i - 18) = dtb(r2, i);
    end
    for i = 23:26
        fuv1(i - 22) = dtb(r1, i);
        fuv2(i - 22) = dtb(r2, i);
    end
    xTy1 = reshape(xTy1, [3, 3]);
    xTy2 = reshape(xTy2, [3, 3]);
    xTz1 = reshape(xTz1, [3, 3]);
    xTz2 = reshape(xTz2, [3, 3]);
    p1 = reshape(p1, [2, 2]);
    p2 = reshape(p2, [2, 2]);
    fuv1 = reshape(fuv1, [2, 2]);
    fuv2 = reshape(fuv2, [2, 2]);
    xTy1 = xTy1';
    xTy2 = xTy2';
    xTz1 = xTz1';
    xTz2 = xTz2';
    p1 = p1';
    p2 = p2';
    fuv1 = fuv1';
    fuv2 = fuv2';
    disp('xTy1')
    disp(xTy1)
    disp('xTz1')
    disp(xTz1)
    disp('p1')
    disp(p1)
    disp('fuv1')
    disp(fuv1)
    disp('xTy2')
    disp(xTy2)
    disp('xTz2')
    disp(xTz2)
    disp('p2')
    disp(p2)
    disp('fuv2')
    disp(fuv2)
    hessi_for_one_letter(p1, fuv1, xTy1, xTz1, lambda, alpha);
    hessi_for_one_letter(p2, fuv2, xTy2, xTz2, lambda, alpha);
    A = zeros(16, 16);
    B = zeros(16, 16);
    hess = zeros(16, 16);
    [A, B, hess] = hessi(p1, p2, fuv1, fuv2, xTy1, xTy2, xTz1, xTz2, lambda, alpha);
    disp('A')
    disp(A)
    disp('B')
    disp(B)
    disp('hessian')
    disp(hess)
    for i = 1:16
        if hess(i, i) < -1e-4
            disp('May not be PSD')
            break;
        end
    end

    fclose(input_file);

end
