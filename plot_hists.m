clear
clc

Lmax_range  = [1 5 10 15];
snr=5;
noise_variance = sqrt(10^(-snr/10));    

Nt=32;Lt=8; Nr=8;K=10;

totalMCrealizations =  1000;
numOfRFchains = zeros(totalMCrealizations, length(Lmax_range));

for k_index = 1:length(Lmax_range)
    Lmax = Lmax_range(k_index);
    for r=1:totalMCrealizations
        [~, ~, numOfRFchains_tmp] = systemModel(Nt, Nr, Lt, K, noise_variance,  Lmax);
        numOfRFchains(r, k_index) = numOfRFchains_tmp;
    end
end

figure;
subplot(1,4,1)
histogram(numOfRFchains(:, 1))
grid on;
xlabel('Number of RF chains', 'FontSize', 11)
ylabel('Frequency', 'FontSize', 11)
title('Lmax=1')
subplot(1,4,2)
histogram(numOfRFchains(:, 2))
grid on;
xlabel('Number of RF chains', 'FontSize', 11)
ylabel('Frequency', 'FontSize', 11)
title('Lmax=5')
subplot(1,4,3)
histogram(numOfRFchains(:, 3))
grid on;
xlabel('Number of RF chains', 'FontSize', 11)
ylabel('Frequency', 'FontSize', 11)
title('Lmax=10')
subplot(1,4,4)
histogram(numOfRFchains(:, 4))
grid on;
xlabel('Number of RF chains', 'FontSize', 11)
ylabel('Frequency', 'FontSize', 11)
title('Lmax=15')

savefig(['./results/hists_K_',num2str(K),'.fig'])
saveas(gcf,['./results/hists_K_',num2str(K),'.eps'],'epsc')