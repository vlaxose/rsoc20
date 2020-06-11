function [rate, power] = systemModel(Nt, Nr, Lt, K, noise_variance)

    %% Initializations
    Prf = 1;
    Pfix = 1;
    rho = 1;
    sigma_phi = 0.01;
    
    %% mmWave channel model
    maxDistance = 5;
    nodeLocations = zeros(K, 2);
    nodeLocations(:, 1) = randi(maxDistance, K, 1);
    nodeLocations(:, 2) = randi(maxDistance, K, 1);
    H = mmWaveChannel(Nr, Nt, K, nodeLocations);
    U = zeros(Nr, Nr, K);
    DD = zeros(Nr, Nt, K);
    V = zeros(Nt, Nt, K);
    Dideal = 1/sqrt(Nt)*fft(eye(Nt));
    Phifull = 1/sqrt(Nt)*sqrt(sigma_phi/2)*(randn(Nt) + 1j*randn(Nt));
    D = Dideal + Phifull;
    for k=1:K
            [U(:,:,k),DD(:,:,k),V(:,:,k)] = svd(H(:,:,k));
            [F_BB(:,:,k),F_RF_ideal(:,:,k)] = hybrid_precoder(V(:,:,k), Nt, Lt, Lt);
            Phi(:,:,k) = Phifull(:,1:Lt);
            F_RF(:,:,k) = F_RF_ideal(:,:,k) + Phi(:,:,k); 
    end
    X = fft(eye(Lt));
    x = X(:,1);
    
    %% Digital beamforming at TX and RX, no interference at the users and beamforming noise
    rate(:, 1) = 0;
    for k=1:K
      eta_k = 0;
      for l=[1:k-1 k+1:K]
        eta_k = eta_k + norm(rho*U(:, 1, l)'*H(:,:,l)*V(:, 1, l))^2;
      end
      rate(:, 1) = rate(:, 1) + real(log2(1+1/(noise_variance^2+eta_k)*norm(rho*U(:, 1, k)'*H(:,:,k)*V(:, :, k))^2));
    end
    power(:, 1) = Pfix + Nt*Prf;
 
    %% Analog beamforming at TX and digital at RX, no interference at the users
    rate(:, 2) = 0;
    for k=1:K
      eta_k = 0;
      for l=[1:k-1 k+1:K]
        eta_k = eta_k + norm(U(:, 1, k)'*H(:,:,l)*D(:, 1))^2;
      end
      rate(:, 2) = rate(:, 2) + real(log2(1+1/(noise_variance^2+eta_k)*norm(U(:, 1, k)'*H(:,:,k)*D(:, 1))^2));    
    end
    power(:, 2) = Pfix + 1*Prf;
    
    %% Hybrid beamforming at TX and digital at RX
    rate(:, 3) = 0;
    for k=1:K
      eta_k = 0;
      for l=[1:k-1 k+1:K]
        eta_k = eta_k + norm(U(:, 1, l)'*H(:,:,l)*F_RF(:,:,l))^2;
      end
      zeta_k = norm(U(:, 1, k)'*H(:,:,k)*Phi(:,:,k))^2;
      rate(:, 3) = rate(:, 3) + real(log2(1+1/(noise_variance^2+eta_k+zeta_k)*norm(U(:, 1, k)'*H(:,:,k)*F_RF_ideal(:,:,k))^2));
      
    end
    power(:, 3) = Pfix + Lt*Prf;
    %% Hybrid beamforming at TX with RF selection and digital at RX, no interference at the users
    % Iterative RF minimization
    ee = 0;
    S_rf_min = zeros(Lt);
    for lt=1:Lt
        S_rf_min(1:lt, 1:lt) = eye(lt);
        rate_rf_min = 0;
        for k=1:K
          eta_k = 0;
          for l=[1:k-1 k+1:K]
            eta_k = eta_k + norm(U(:, 1, l)'*H(:,:,l)*D(:, 1:Lt)*S_rf_min)^2;
          end
          zeta_k = norm(U(:, 1, k)'*H(:,:,k)*Phi(:,:,k)*S_rf_min)^2;
          rate_rf_min = rate_rf_min + real(log2(1+1/(noise_variance^2+eta_k+zeta_k)*norm(U(:, 1, k)'*H(:,:,k)*Dideal(:, 1:Lt)*S_rf_min)^2));
        end
        power_rf_min = Pfix + lt*Prf;
        ee_rf_min = rate_rf_min/power_rf_min;
        if ee_rf_min > ee
            rate(:, 4) = rate_rf_min; 
            power(:, 4) = power_rf_min;
            ee = ee_rf_min;
        end
    end
    % Exhaustive search best RF subset selection
    ee = 0;
    S_bf_opt =[];
    if(Lt<=12)
        binary_combinations = dec2bin(0:2^Lt - 1) - '0';        
        for lt=1:size(binary_combinations, 1)
           S_bf = diag(binary_combinations(lt, :));
           rate_rf_bf = 0;
           for k=1:K
             eta_k = 0;
             for l=[1:k-1 k+1:K]
               eta_k = eta_k + norm(U(:, 1, l)'*H(:,:,l)*D(:,1:Lt)*S_bf)^2;
             end
             zeta_k = norm(U(:, 1, k)'*H(:,:,k)*Phi(:,:,k)*S_bf)^2;
             rate_rf_bf = rate_rf_bf +  real(log2(1+1/(noise_variance^2+eta_k+zeta_k)*norm(U(:, 1, k)'*H(:,:,k)*Dideal(:,1:Lt)*S_bf)^2));
           end
           power_rf_bf = Pfix + sum(diag(S_bf))*Prf;
           ee_rf_bf = rate_rf_bf/power_rf_bf;
           if ee_rf_bf>ee
               S_bf_opt = S_bf;
               rate(:, 5) = rate_rf_bf;
               power(:, 5) = power_rf_bf;
               ee = ee_rf_bf;
           end
        end
    else
        rate(:, 5) = 0;
        power(:, 5) = 0;
    end
diag(S_bf_opt)'

    Pe = 0.1*Pfix;
    S = proposed(H, U, V, D, Dideal, Phifull, K, Nt, Lt, noise_variance, Prf, Pfix+Pe);
diag(S)'
    rate(:, 6) = 0;
    for k=1:K
      eta_k = 0;
      for l=[1:k-1 k+1:K]
        eta_k = eta_k + norm(U(:,1,l)'*H(:,:,l)*D*S)^2;
      end
      zeta_k = norm(U(:, 1, k)'*H(:,:,k)*Phifull*S)^2;
      rate(:, 6) = rate(:, 6) + real(log2(1+1/(noise_variance^2+eta_k+zeta_k)*norm(U(:,1,k)'*H(:,:,k)*Dideal*S)^2));
      
    end
    power(:, 6) = (Pfix+Pe) + sum(diag(S))*Prf;
    
    
end
