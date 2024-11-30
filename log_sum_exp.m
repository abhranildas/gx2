function s=log_sum_exp(logs,dim)
% log sum exp of a batch of logs (base 10), where the smallest is no less than
% realmin=1e-308 of the largest
if nargin <2
    dim='all';
end

if isempty(logs)
    s=-inf;
else
    % logs=sort(logs,dim,'descend');
    m=max(logs,[],dim);
    % m=logs(1); % largest log
    s=m+log10(sum(10.^(logs-m),dim));
end
end