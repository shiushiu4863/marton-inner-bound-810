function result = random_channel(row, column)
    result = [];
    for i = 1:row
        result = vertcat(result, random_vector(column));
    end
end