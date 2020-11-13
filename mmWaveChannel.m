function H = mmWaveChannel(Nr, Nt, K, nodeLocations, Lmax)

% Initialization
H = zeros(Nr, Nt, K);
Ghz = 28;
wavelength = 30/Ghz; % w=c/lambda

theta = zeros(K,1);
phi = zeros(K,1);
distanceFromBS = zeros(K, 2);
for k=1:K
  theta(k) = atan(nodeLocations(k, 2, 1)/nodeLocations(k, 1, 1));
  phi(k) = pi/2 - theta(k);
  distanceFromBS(k, 1) = sqrt(nodeLocations(k, 2, 1)^2 + nodeLocations(k, 1, 1)^2);
end

for k=1:K
  Lk = randi(Lmax);
	for ell=1:Lk
    ar = angle(theta(k)*randn, Nr, wavelength);
   	at = angle(phi(k)*randn, Nt, wavelength);
 		if(ell==1)
			channel_coef = 1*randn(1)/(4*pi/wavelength*distanceFromBS(k, 1)^2);
		else
			channel_coef = 0.8*randn(1)/(4*pi/wavelength*distanceFromBS(k, 1)^2);

    end
		H(:,:,k) = H(:,:,k) + channel_coef * ar*at';
  end
	H(:,:,k) = sqrt(Nr*Nt)/sqrt(Lk)*H(:,:,k);
end
end

% Generate the transmit and receive array responces
function vectors_of_angles=angle(phi, M, wavelength)

    array_element_spacing = 0.5*wavelength;
    wavenumber = 2*pi/wavelength; % k = 2pi/lambda
    phi0 = 0; % mean AOA
    phase_shift = wavenumber*array_element_spacing*sin(phi0-phi)*(0:M-1).';
    vectors_of_angles = 1/sqrt(M)*exp(1j*phase_shift);
end
