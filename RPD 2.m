function vHat = decodeLogDomain_RPD_correct(rx, H, N0, max_iterations)
    [M, N] = size(H);
    C_v = (4*rx./N0);
    
    R_mn = zeros(M, N);      % Current iteration
    R_mn_prev = zeros(M, N); % Previous iteration ← NEED THIS
    Q_n = C_v;
    Q_n_prev = C_v;          % Previous posterior ← NEED THIS
    
    for iter = 1:max_iterations
        
        % Use PREVIOUS iteration for decisions
        x_hat = double(Q_n_prev < 0);
        S = mod(x_hat' * H', 2);
        
        if all(S == 0)
            vHat = x_hat';
            return;
        end
        
        % Compute RP
        f = zeros(1, N);
        for n = 1:N
            M_n = find(H(:, n));
            f(n) = sum(S(M_n));
        end
        
        [~, order] = sort(f, 'descend');
        N_star = sum(f > 0);
        
        if N_star == 0
            vHat = x_hat';
            return;
        end
        
        % ===== EQUATION (6): Initialize current iteration =====
        Q_n = Q_n_prev;  % Q_n^(i) = Q_n^(i-1)
        
        % Store previous iteration messages
        R_mn_prev = R_mn;  % Save R^(i-1)
        
        % Update unreliable nodes
        for idx = 1:N_star
            n = order(idx);
            M_n = find(H(:, n));
            M_n_sorted = sort(M_n); % Sort for proper ordering
            
            for m = M_n_sorted'
                
                % ===== EQUATION (7): V2C message =====
                Q_mn = Q_n(n);
                
                for m_prime = M_n_sorted'
                    if m_prime ~= m
                        if m_prime < m
                            % Already updated in CURRENT iteration
                            Q_mn = Q_mn - R_mn(m_prime, n);
                        else
                            % Not yet updated, use PREVIOUS iteration
                            Q_mn = Q_mn - R_mn_prev(m_prime, n);
                        end
                    end
                end
                
                % ===== EQUATION (8): C2V message =====
                N_m = find(H(m, :));
                product = 1;
                
                for n_prime = N_m
                    if n_prime ~= n
                        M_n_prime = find(H(:, n_prime));
                        M_n_prime_sorted = sort(M_n_prime);
                        Q_mn_prime = Q_n(n_prime);
                        
                        for m_pp = M_n_prime_sorted'
                            if m_pp ~= m
                                if m_pp < m
                                    Q_mn_prime = Q_mn_prime - R_mn(m_pp, n_prime);
                                else
                                    Q_mn_prime = Q_mn_prime - R_mn_prev(m_pp, n_prime);
                                end
                            end
                        end
                        product = product * tanh(Q_mn_prime / 2);
                    end
                end
                
                product = max(min(product, 1-1e-10), -1+1e-10);
                R_mn(m, n) = 2 * atanh(product);
                
                % ===== EQUATION (9): Update posterior =====
                Q_n(n) = C_v(n) + sum(R_mn(M_n_sorted, n));
            end
        end
        
        % Save for next iteration
        Q_n_prev = Q_n;
    end
    
    vHat = double(Q_n < 0)';
end