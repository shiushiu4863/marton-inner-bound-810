function result = random_vector(len)
    result = zeros(1, len);
    for i = 1:len - 1
        result(i) = (1 - sum(result)) * rand;
    end
    result(len) = 1 - sum(result);
end