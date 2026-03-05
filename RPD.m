function vHat = decodeLogDomain_RPD(rx, H, N0, max_iterations)
    [M, N] = size(H);
    
    % Prior log-likelihood (Equation 5)
    % For BPSK: LLR = 4*rx/N0 (we received rx, sent tx where tx = 1-2*bit)
    C_v = (4*rx./N0);  % Note: positive LLR means bit = 0
    
    % Initialize messages
    R_mn = zeros(M, N);     % Check to variable
    Q_n = C_v;              % Posterior LLR
    
    % Main iteration loop
    for iter = 1:max_iterations
        
        % ===== COMPUTE RELIABILITY PROFILE =====
        % Hard decision (LLR > 0 means bit = 0)
        x_hat = double(Q_n < 0);  % Note: < 0 means bit = 1
        
        % Syndrome S = x_hat * H^T (mod 2)
        S = mod(x_hat' * H', 2);
        
        % Check convergence
        if all(S == 0)
            vHat = x_hat';
            return;
        end
        
        % Compute reliability for each variable node
        f = zeros(1, N);
        for n = 1:N
            M_n = find(H(:, n));
            f(n) = sum(S(M_n));
        end
        
        % Sort by reliability (descending)
        [~, order] = sort(f, 'descend');
        N_star = sum(f > 0);
        
        if N_star == 0
            vHat = x_hat';
            return;
        end
        
        % ===== UPDATE UNRELIABLE NODES IN ORDER =====
        for idx = 1:N_star
            n = order(idx);
            
            M_n = find(H(:, n));
            
            for m = M_n'
                % Compute V2C message (Equation 7)
                Q_mn = Q_n(n);
                for m_prime = M_n'
                    if m_prime ~= m
                        Q_mn = Q_mn - R_mn(m_prime, n);
                    end
                end
                
                % Update C2V message (Equation 8)
                N_m = find(H(m, :));
                product = 1;
                
                for n_prime = N_m
                    if n_prime ~= n
                        % Compute V2C for n_prime
                        M_n_prime = find(H(:, n_prime));
                        Q_mn_prime = Q_n(n_prime);
                        for m_pp = M_n_prime'
                            if m_pp ~= m
                                Q_mn_prime = Q_mn_prime - R_mn(m_pp, n_prime);
                            end
                        end
                        product = product * tanh(Q_mn_prime / 2);
                    end
                end
                
                % Numerical stability
                product = max(min(product, 1-1e-10), -1+1e-10);
                R_mn(m, n) = 2 * atanh(product);
                
                % Update posterior (Equation 9)
                Q_n(n) = C_v(n) + sum(R_mn(M_n, n));
            end
        end
    end
    
    % Final decision
    vHat = double(Q_n < 0)';
end