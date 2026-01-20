% function l2error = L2Error(a, b)

% if ~isequal(size(a), size(b))
%     error('Vectors must have the same size.');
% end

% N = numel(a);
% sum_sq = sum((a - b).^2);
% l2error = sqrt(sum_sq / N);

% end

function l2error = L2Error(a, b)
    if ~isequal(size(a), size(b))
        error('Vectors/Matrices must have the same size.');
    end
    
    N = numel(a);
    sum_sq = sum((a(:) - b(:)).^2); 
    
    l2error = sqrt(sum_sq / N);
end