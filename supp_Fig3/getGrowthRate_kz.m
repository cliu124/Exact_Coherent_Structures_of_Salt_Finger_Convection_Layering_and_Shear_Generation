function GR = getGrowthRate_kz(kz_list,Ra,Pr,Rrho,tau,aw0)
    % This function uses Parallel Computing Toolbox to rise the computing
    % speed. Normally, a single core will be used when Parallel Computing 
    % Toolbox is not installed on the local machine.
    sizeOfkz = size(kz_list,1);
    % create an array to store values of growth rate
    GR = zeros(sizeOfkz,1);
    if canUseParallelPool
        % Parallel Computing Toolbox is installed
        parfor kz_index=1:sizeOfkz
            kz=kz_list(kz_index);
            % compute eigenvalues
            [~,eigenvalue] = eig_Nguyen2026(aw0,Ra,Pr,Rrho,tau,kz);
            % eig_val(find(abs(eig_val)>10^1))=-Inf;
            clean_mask = ~isinf(eigenvalue) & ~isnan(eigenvalue) & (abs(eigenvalue) < 1e2);
            eigenvalue = eigenvalue(clean_mask);
            [~, sort_indices] = sort(real(eigenvalue), 'descend');
            eigenvalue = eigenvalue(sort_indices);
            % compute growth rate
            GR(kz_index)=eigenvalue(1);
        end
    else
        % Parallel Computing Toolbox is not installed
        for kz_index=1:sizeOfkz
            kz=kz_list(kz_index);
            % compute eigenvalues
            [~,eigenvalue] = eig_Nguyen2026(aw0,Ra,Pr,Rrho,tau,kz);
            % eig_val(find(abs(eig_val)>10^1))=-Inf;
            clean_mask = ~isinf(eigenvalue) & ~isnan(eigenvalue) & (abs(eigenvalue) < 1e2);
            eigenvalue = eigenvalue(clean_mask);
            [~, sort_indices] = sort(real(eigenvalue), 'descend');
            eigenvalue = eigenvalue(sort_indices);
            % compute growth rate
            GR(kz_index)=eigenvalue(1);
        end
    end
end