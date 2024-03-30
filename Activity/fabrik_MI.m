function [p] = fabrik_MI(p, n, t, tol)
    d = sqrt(sum(diff(p, 1).^2, 2));
    
    dist = norm(p(1,:) - t);
    
    if dist > sum(d)
        disp('Target is unreachable!')
        lambdas = d ./ sqrt(sum((t - p(1:end-1,:)).^2, 2));
        p(2:end,:) = (1 - lambdas) .* p(1:end-1,:) + lambdas .* t;
    
    else
        disp('Target is reachable!')
        b = p(1,:);
        difa = norm(p(n,:) - t);
    
        while difa > tol
            p(n,:) = t;
            lambdas = d ./ sqrt(sum(diff(p(n-1:-1:1,:), 1).^2, 2));
            p(1:end-1,:) = (1 - lambdas) .* p(n-1:-1:1,:) + lambdas .* p(1:end-1,:);
    
            p(1,:) = b;
            lambdas = d ./ sqrt(sum(diff(p(1:end-1,:), 1).^2, 2));
            p(2:end,:) = (1 - lambdas) .* p(1:end-1,:) + lambdas .* p(2:end,:);
            difa = norm(p(n,:) - t);
        end
    end
    p
end