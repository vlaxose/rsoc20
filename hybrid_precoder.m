function [F_BB,F_RF] = hybrid_precoder(V, Nt, Lt, Ns)

    F_BB = zeros(Lt, Ns);
    F_RF = zeros(Nt, Lt);
    D = 1/sqrt(Nt)*fft(eye(Nt));
    for indx=1:100
      % Get the index with the maximum energy
      diff = V(:,1:Ns) - F_RF*F_BB;
      Psi = D'*diff/norm(diff, 'fro');
      C = diag(Psi*Psi');
      [~, I] = sort(C, 'descend');

      % Update the precoders
      F_t = [D(:,I(1)) F_RF];
      F_RF = F_t(:,1:Lt);
      F_BB = pinv(F_RF)*V(:,1:Ns);

    end
    F = F_RF*F_BB;
%     F_BB = rho*F_BB/norm(F, 'fro');

   
    if(norm(V(:,1:Ns) - F)^2/norm(V)^2>1e-1)
        warning(['A/D beamformer did not converge, with error:', num2str(norm(V(:,1:Ns) - F, 'fro'))]);
    end

end