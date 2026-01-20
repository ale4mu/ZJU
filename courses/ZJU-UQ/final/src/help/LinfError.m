function linf_error = LinfError(a, b)

if ~isequal(size(a), size(b))
    error('Vectors must have the same size.');
end

max_diff = max(abs(a - b));
linf_error = max_diff;

end