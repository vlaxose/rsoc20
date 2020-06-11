function [ee, rate, power] = ee_computation(S, noise_variance, Nr, H, F_RF, F_BB, W, Prf)


    % Noise covariance matrix
    K = real(W'*H*F_RF*S*(F_BB*F_BB')*S*F_RF'*H'*W + noise_variance^2*eye(Nr));
    
    % Achievable information rate
    rate = real(log2(det(eye(Nr) + 1/Nr*pinv(K)*W'*H*F_RF*S*(F_BB*F_BB')*S'*F_RF'*H'*W)));
    
    % Consumed power
    F = F_RF*S*F_BB;
    power = trace(S)*Prf;
    
    % Energy efficiency computation
    ee = rate/power;
end
