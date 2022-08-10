function result = h(p)
    result = - sum(p .* log2(p), 'all', 'omitnan');
end