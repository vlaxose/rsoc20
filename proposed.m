function S = proposed(H, U, V, D, Dideal, Phifull, K, Nt, Lt, noise_variance, Prf, Pfix)

  maxDBiters = 5;

% 
%   Rmin = 0;
% %   R = 0;
%   for k=1:K
% % % % 
% % % %     eta_k = 0;
% % % %     for l=[1:k-1 k+1:K]
% % % %       eta_k = eta_k + norm(U(:, 1, l)'*H(:,:,l)*F_RF*F_BB(:, 1))^2;
% % % %     end
%       Rmin = Rmin + real(log2(1+1/(noise_variance^2)*norm(U(:, 1, k)'*H(:,:,k)*V(:, 1, k))^2));
% %         R = R + real(trace(1/(noise_variance^2)* rho^2*U(:, 1, k)'*H(:,:,k)*H(:,:,k)'*U(:, 1, k)));
% %       
%   end
% %   Rmin
% %   R
  
  kappa = zeros(1, maxDBiters+1);
  indx = 1;
  % Dinkelbach iterations
  for iter=1:maxDBiters 

    % CVX solution
    cvx_begin quiet
      variable s(Nt)
      R = 0;
      for k=1:K
        eta_k = 0;
        for l=[1:k-1 k+1:K]
          eta_k = eta_k + real(trace(U(:, 1, l)'*H(:,:,l)*D*diag(s)*D'*H(:,:,l)'*U(:, 1, l)));
        end
        if(norm(Phifull)>0)
          zeta_k = real(trace(U(:, 1, k)'*H(:,:,k)*Phifull*diag(s)*Phifull'*H(:,:,k)'*U(:,1,k)));
        else
          zeta_k = 0;
        end
        R = R + real(trace(1/noise_variance^2*U(:, 1, k)'*H(:,:,k)*Dideal*diag(s)*Dideal'*H(:,:,k)'*U(:, 1, k)));
        

      end

      P = Pfix+sum(s)*Prf;
                

      maximize(R-(kappa(iter)*P+eta_k+zeta_k))

      subject to
         s >= 0;
         s <= 1;
%          R >= 0.8*Rmin;
         P <= Prf*Lt + Pfix;
    cvx_end

    % Thresholding to transform from real to binary
    s_indx = find(s/norm(s)>1e-6);
    Lt_opt = length(s_indx);
    S = zeros(Nt);
    S(s_indx(1:Lt_opt), s_indx(1:Lt_opt)) = eye(Lt_opt);

    % Check the solution
    if(isnan(cvx_optval) || isinf(cvx_optval) || sum(s<0) == length(s) || length(find(diag(S)>0)) == length(s))
      S = zeros(Nt);
      S(1:1, 1:1) = eye(1);
    end
    
    for k=1:K
      eta_k = 0;
      for l=[1:k-1 k+1:K]
          eta_k = eta_k + real(trace(U(:, 1, l)'*H(:,:,l)*D*S*D'*H(:,:,l)'*U(:, 1, l)));
      end
      zeta_k = norm(U(:, 1, k)'*H(:,:,k)*Phifull*S)^2;

%         R = R +  real(log2(1+1/(noise_variance^2+eta_k)*norm(U(:, 1, k)'*H(:,:,k)*D*S*x)^2));
        R = R + real(trace(1/(noise_variance^2+eta_k+zeta_k)*U(:, 1, k)'*H(:,:,k)*Dideal*S*Dideal'*H(:,:,k)'*U(:, 1, k)));

      
    end
    kappa(indx+1) =  R/(Pfix+trace(S)*Prf);

    
    indx = indx + 1;

    
  end



end
