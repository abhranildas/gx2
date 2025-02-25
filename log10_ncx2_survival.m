function log10SF = log10_ncx2_survival(x, nu, lambda)
% log10_ncx2_survival 
%     Compute log10 of the survival function (1 - CDF) of the noncentral
%     chi-square distribution with degrees of freedom nu and
%     noncentrality lambda at point x.
%
%     That is, returns log10( Q_nc(x; nu, lambda) ), where
%     Q_nc = 1 - F_{ncx^2}(x; nu, lambda).
%
% Usage:
%     log10SF = log10_ncx2_survival(x, nu, lambda)
%
% Inputs:
%     x        - Scalar or vector of evaluation points (x>0)
%     nu       - Degrees of freedom (nu>0)
%     lambda   - Noncentrality parameter (lambda>=0)
%
% Output:
%     log10SF  - log10( Survival function ) of the noncentral chi-square
%                distribution at x.

    % --- Input checks (optional) ---
    if any(x <= 0)
        error('All x must be > 0.');
    end
    if nu <= 0
        error('nu must be > 0.');
    end
    if lambda < 0
        error('lambda must be >= 0.');
    end

    % Reshape x into a column vector for convenience
    x = x(:);
    log10SF = zeros(size(x));

    % If lambda == 0, this is the *central* chi-square distribution,
    % so the survival function is simply gammainc(x/2, nu/2, 'upper').
    % We handle that case separately to avoid log10(0) issues.
    if lambda == 0
        % For each x(i), compute Q_upper = gammainc(x(i)/2, nu/2, 'upper'),
        % then take log10 of it. If Q_upper is 0, we set -Inf.
        for i = 1:numel(x)
            val = gammainc(x(i)/2, nu/2, 'upper');  % regularized inc. gamma
            if val <= 0
                log10SF(i) = -Inf;  % Probability is numerically zero
            else
                log10SF(i) = log10(val);
            end
        end
        
        return;  % We are done with lambda=0 case
    end

    % --- Otherwise, use the noncentral series expansion ---
    
    % Series parameters
    maxIter = 2000;      % max # of terms
    reltol  = 1e-14;     % relative tolerance

    for i = 1:numel(x)
        xi = x(i);

        % Initialize partial sum in log10 domain: log10(0) = -Inf
        partialSum_log10 = -Inf;

        for k = 0:maxIter
            % ---------------------------
            % 1) Compute log10 of p_k:
            %
            %    p_k = exp(-lambda/2) * (lambda/2)^k / k!
            %
            %    log10(p_k) = - (lambda/2)*1/log(10) + k*log10(lambda/2) - log10(k!)
            %               = - (lambda/2)/ln(10) + k*log10(lambda/2) - (gammaln(k+1)/ln(10))
            % ---------------------------
            log10_p_k = - (lambda/2) / log(10) ...
                        + k * log10(lambda/2) ...
                        - (gammaln(k+1) / log(10));

            % ---------------------------
            % 2) Q_k = gammainc(xi/2, (nu+2*k)/2, 'upper')
            %    log10_Q_k = log10(Q_k). If Q_k = 0 => -Inf
            % ---------------------------
            Q_k = gammainc(xi/2, (nu + 2*k)/2, 'upper');
            if Q_k <= 0
                log10_Q_k = -Inf;
            else
                log10_Q_k = log10(Q_k);
            end

            % ---------------------------
            % 3) Term = p_k * Q_k => log10_term = log10_p_k + log10_Q_k
            % ---------------------------
            log10_term = log10_p_k + log10_Q_k;

            % ---------------------------
            % 4) Stable sum in base 10
            % ---------------------------
            new_partialSum_log10 = logplus10(partialSum_log10, log10_term);

            % ---------------------------
            % 5) Check for convergence
            %    Relative change ~ 10^(new - old) - 1
            % ---------------------------
            if new_partialSum_log10 > partialSum_log10
                delta = 10^(new_partialSum_log10 - partialSum_log10) - 1;
            else
                delta = 0;
            end

            partialSum_log10 = new_partialSum_log10;

            if delta < reltol
                break;  % Series has converged enough
            end
        end

        % partialSum_log10 is log10(SF) for x(i)
        log10SF(i) = partialSum_log10;
    end
end

function s = logplus10(a, b)
% logplus10(a, b) = log10(10^a + 10^b) in a stable way.
%
%   If a = -Inf or b = -Inf, just return the other value.
%   Otherwise, do base-10 "log-sum-exp" trick.

    if isinf(a) && a < 0
        s = b;
    elseif isinf(b) && b < 0
        s = a;
    elseif a > b
        s = a + log10(1 + 10^(b - a));
    else
        s = b + log10(1 + 10^(a - b));
    end
end
