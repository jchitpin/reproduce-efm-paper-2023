function S = simulate_stoich(n)
    
    % Unique reaction pairs (precursor to product)
    i = 0;
    rxn_pairs = zeros(n(2), 2);
    while i < n(2)
        rxn_rand = datasample(1:n(1), 2, 'Replace', false); % no replacement to avoid same index for precursor/product
        if ~ismember(rxn_rand, rxn_pairs, 'rows')
            i = i + 1;
            rxn_pairs(i,1) = rxn_rand(1);
            rxn_pairs(i,2) = rxn_rand(2);
        end
    end
    rxn_pairs(:,1) = rxn_pairs(:,1) * -1;
    
    % Check that each metabolite index is a precursor and product
    precursor_needed = find(~ismember([-(1:n(1))], rxn_pairs(:,1)));
    product_needed = find(~ismember([1:n(1)], rxn_pairs(:,2)));
    
    % If there are empty metabolite indices, add them back in
    if ~isempty(precursor_needed) || ~isempty(product_needed)
        
        num_prune = max(length(precursor_needed), length(product_needed));
    
        % Prune reaction pairs with extra precursors/products
        [a1,b1] = hist(rxn_pairs(:,1), unique(rxn_pairs(:,1)));
        [a2,b2] = hist(rxn_pairs(:,2), unique(rxn_pairs(:,2)));
        rm_precursor = b1(a1 > 1);
        rm_product = b2(a2 > 1);

        intersect1 = ismember(rxn_pairs(:,1), rm_precursor);
        intersect2 = ismember(rxn_pairs(:,2), rm_product);

        combined = intersect1 + intersect2;
        combined_idx = find(combined > 1); % indices of rxns that can be removed

        if length(combined_idx) < num_prune
            error('Error. Try running this function again with a different random seed.')
        end

        % Must ensure pairs with the same precursors/products are not exhausted
        remove_idx = [];
        remove_pairs = [0 0];
        while length(remove_idx) < num_prune
            temp_idx = datasample(combined_idx, 1);
            %combined_idx = combined_idx(setdiff(1:end, temp_idx));
            temp_pair = rxn_pairs(temp_idx,:);
            if ~ismember(temp_pair, remove_pairs, 'rows')
                remove_pairs = vertcat(remove_pairs, temp_pair);
                remove_idx = vertcat(remove_idx, temp_idx);
            end
        end
        rxn_pairs = rxn_pairs(setdiff(1:end, remove_idx), :);

        % Add reaction pairs with missing precursors/products
        % Must ensure the same index is both precursor/product
        num_add = min(length(precursor_needed), length(product_needed));
        precursor = precursor_needed(1:num_add);
        product = product_needed(1:num_add);
        brk = 0;
        while any(precursor == product)
            brk = brk + 1;
            precursor = precursor(randperm(length(precursor)));
            product = product(randperm(length(product)));
            if brk == 500
               error('Error. Try again with a different seed.') 
            end
        end

        %rxn_pairs = [rxn_pairs; [-precursor_needed(1:num_add); product_needed(1:num_add)]'];
        rxn_pairs = [rxn_pairs; [-precursor; product]'];
        
        if num_prune > num_add
           % Update which metabolite indices are missing
            precursor_needed = find(~ismember([-(1:n(1))], rxn_pairs(:,1)));
            product_needed = find(~ismember([1:n(1)], rxn_pairs(:,2)));
            
            if ~isempty(precursor_needed)
                rxn_pairs = [rxn_pairs; [-precursor_needed; datasample(setdiff(1:n(1), precursor_needed), length(precursor_needed), 'Replace', true)]'];          
            end
            if ~isempty(product_needed)
                rxn_pairs = [[-datasample(setdiff(1:n(1), product_needed), length(product_needed), 'Replace', true); product_needed]'; rxn_pairs];          
            end
        end
       
        % Error-checking
        if size(rxn_pairs,1) ~= n(2)
            error('Error. Number of reactions (%s) does not match n(2).', num2str(size(rxn_pairs,1)))
        end
        if isempty(~ismember([-(1:n(1))], rxn_pairs(:,1)))
            error('Error. Network is not fully connected.')
        end
        if isempty(~ismember([1:n(1)], rxn_pairs(:,2)))
            error('Error. Network is not fully connected.')
        end
    end
    
    % Convert reaction pairs to stoichiometry matrix
    S = zeros(n(1), n(2));
    for i = 1:n(2)
        S(-rxn_pairs(i,1), i) = -1;
        S(+rxn_pairs(i,2), i) = +1;
    end
    
    if ~all(size(unique(S', 'rows')) == size(S'))
        error('Error. Duplicate reactions in stoichiometry matrix not allowed.')
    end
    if any(sum(S, 1) ~= 0)
        change_idx = find((rxn_pairs(:,1) + rxn_pairs(:,2)) == 0);
        if change_idx < 1
            change_idx(change_idx < 1) = 2;
        end
        rxn_pairs(change_idx,2) = rxn_pairs(change_idx,2) - 1;
        S = zeros(n(1), n(2));
        for i = 1:n(2)
            S(-rxn_pairs(i,1), i) = -1;
            S(+rxn_pairs(i,2), i) = +1;
        end
    end
    if any(sum(S, 1) ~= 0)
        change_idx = find((rxn_pairs(:,1) + rxn_pairs(:,2)) == 0);
        if change_idx < 1
            change_idx(change_idx < 1) = 2;
        end
        rxn_pairs(change_idx,2) = rxn_pairs(change_idx,2) - 1;
        disp(rxn_pairs)
        error('Error. Each reaction must contain one precursor and one product.')
    end
    S
end
