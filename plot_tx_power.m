clear
clc

% snr_range = [0 3 5 8 11 14 17 20];
snr_range = [0 5 11 17];
Nt=8;Lt=8;Nr=8;K=5;Lmax=15;

totalMCrealizations = 100;
rate = zeros(totalMCrealizations, length(snr_range), 6);
power = zeros(totalMCrealizations, length(snr_range), 6);

for snr_index = 1:length(snr_range)
    noise_variance = sqrt(10^(-snr_range(snr_index)/10));    
    for r=1:totalMCrealizations
        [rate(r, snr_index, :), power(r, snr_index, :)] = systemModel(Nt, Nr, Lt, K, noise_variance, Lmax);
    end
end
mean_rate = squeeze(sum(rate, 1))/totalMCrealizations;
mean_power = squeeze(sum(power, 1))/totalMCrealizations;    

mean_ee = zeros(size(mean_rate));
mean_ee(:, 1) = mean_rate(:, 1)./mean_power(:, 1);
mean_ee(:, 2) = mean_rate(:, 2)./mean_power(:, 2);
mean_ee(:, 3) = mean_rate(:, 3)./mean_power(:, 3);
mean_ee(:, 4) = mean_rate(:, 4)./mean_power(:, 4);
mean_ee(:, 5) = mean_rate(:, 5)./mean_power(:, 5);
mean_ee(:, 6) = mean_rate(:, 6)./mean_power(:, 6);


figure;
p1=plot(snr_range, mean_rate(:, 1)); hold on;
set(p1,'LineWidth',1, 'LineStyle', '-', 'Color', 'Black');
p2=plot(snr_range, mean_rate(:, 2)); hold on;
set(p2,'LineWidth',1, 'LineStyle', '--', 'Color', 'Black');
p3=plot(snr_range, mean_rate(:, 3)); hold on;
set(p3,'LineWidth',1, 'LineStyle', '-.', 'Color', 'Black');
p4=plot(snr_range, mean_rate(:, 4)); hold on;
set(p4,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Green', 'Marker', 'o');
p5=plot(snr_range, mean_rate(:, 5)); hold on;
set(p5,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Blue', 'Marker', 's');
p6=plot(snr_range, mean_rate(:, 6)); hold on;
set(p6,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Red', 'Marker', 'h');
grid on;
xlabel('SNR (dB)', 'FontSize', 11)
ylabel({'Rate', '(Mbits/sec)'}, 'FontSize', 11)
lg = legend('Digital BF', 'Analog BF', 'Hybrid BF', 'Hybrid BF with RF minimization', 'Hybrid BF with RF exhaustive search', 'Proposed', 'Location', 'Best');
lg.FontSize = 8;

savefig(['./results/rate_snr_with_eta_and_zeta_',num2str(Nt),'_',num2str(Lt),'_',num2str(Nr),'.fig'])
saveas(gcf,['./results/rate_snr_with_eta_and_zeta_',num2str(Nt),'_',num2str(Lt),'_',num2str(Nr),'.eps'],'epsc')

figure;
mean_ee(:, 1) = mean_rate(:, 1)./mean_power(:, 1);
p1=plot(snr_range, mean_ee(:, 1)); hold on;
set(p1,'LineWidth',1, 'LineStyle', '-', 'Color', 'Black');
mean_ee(:, 2) = mean_rate(:, 2)./mean_power(:, 2);
p2=plot(snr_range, mean_ee(:, 2)); hold on;
set(p2,'LineWidth',1, 'LineStyle', '--', 'Color', 'Black');
mean_ee(:, 3) = mean_rate(:, 3)./mean_power(:, 3);
p3=plot(snr_range, mean_ee(:, 3)); hold on;
set(p3,'LineWidth',1, 'LineStyle', '-.', 'Color', 'Black');
mean_ee(:, 4) = mean_rate(:, 4)./mean_power(:, 4);
p4=plot(snr_range, mean_ee(:, 4)); hold on;
set(p4,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Green', 'Marker', 'o');
mean_ee(:, 5) = mean_rate(:, 5)./mean_power(:, 5);
p5=plot(snr_range, mean_ee(:, 5)); hold on;
set(p5,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Blue', 'Marker', 's');
mean_ee(:, 6) = mean_rate(:, 6)./mean_power(:, 6);
p6=plot(snr_range, mean_ee(:, 6)); hold on;
set(p6,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Red', 'Marker', 'h');
grid on;
xlabel('SNR (dB)', 'FontSize', 11)
ylabel({'Energy Efficiency', '(Mbits/Joule)'}, 'FontSize', 11)
lg = legend('Digital BF', 'Analog BF', 'Hybrid BF', 'Hybrid BF with RF minimization', 'Hybrid BF with RF exhaustive search', 'Proposed', 'Location', 'Best');
lg.FontSize = 8;

savefig(['./results/ee_snr_with_eta_and_zeta_',num2str(Nt),'_',num2str(Lt),'_',num2str(Nr),'.fig'])
saveas(gcf,['./results/ee_snr_with_eta_and_zeta_',num2str(Nt),'_',num2str(Lt),'_',num2str(Nr),'.eps'],'epsc')


figure;
p1=plot(snr_range, mean_power(:, 1)); hold on;
set(p1,'LineWidth',1, 'LineStyle', '-', 'Color', 'Black');
p2=plot(snr_range, mean_power(:, 2)); hold on;
set(p2,'LineWidth',1, 'LineStyle', '--', 'Color', 'Black');
p3=plot(snr_range, mean_power(:, 3)); hold on;
set(p3,'LineWidth',1, 'LineStyle', '-.', 'Color', 'Black');
p4=plot(snr_range, mean_power(:, 4)); hold on;
set(p4,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Green', 'Marker', 'o');
p5=plot(snr_range, mean_power(:, 5)); hold on;
set(p5,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Blue', 'Marker', 's');
p6=plot(snr_range, mean_power(:, 6)); hold on;
set(p6,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Red', 'Marker', 'h');
grid on;
xlabel('SNR (dB)', 'FontSize', 11)
ylabel({'Power', '(Joule)'}, 'FontSize', 11)
lg = legend('Digital BF', 'Analog BF', 'Hybrid BF', 'Hybrid BF with RF minimization', 'Hybrid BF with RF exhaustive search', 'Proposed', 'Location', 'Best');
lg.FontSize = 8;

savefig(['./results/power_snr_with_eta_and_zeta_',num2str(Nt),'_',num2str(Lt),'_',num2str(Nr),'.fig'])
saveas(gcf,['./results/power_snr_with_eta_and_zeta_',num2str(Nt),'_',num2str(Lt),'_',num2str(Nr),'.eps'],'epsc')


figure;
p1=plot(mean_rate(:, 1), mean_ee(:, 1)); hold on;
set(p1,'LineWidth',1, 'LineStyle', '-', 'Color', 'Black');
p2=plot(mean_rate(:, 2), mean_ee(:, 2)); hold on;
set(p2,'LineWidth',1, 'LineStyle', '--', 'Color', 'Black');
p3=plot(mean_rate(:, 3), mean_ee(:, 3)); hold on;
set(p3,'LineWidth',1, 'LineStyle', '-.', 'Color', 'Black');
p4=plot(mean_rate(:, 4), mean_ee(:, 4)); hold on;
set(p4,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Green', 'Marker', 'o');
p5=plot(mean_rate(:, 5), mean_ee(:, 5)); hold on;
set(p5,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Blue', 'Marker', 's');
p6=plot(mean_rate(:, 6), mean_ee(:, 6)); hold on;
set(p6,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Red', 'Marker', 'h');
grid on;
xlabel({'Rate', '(Mbits/sec)'}, 'FontSize', 11)
ylabel({'Energy Efficiency', '(Mbits/Joule)'}, 'FontSize', 11)
lg = legend('Digital BF', 'Analog BF', 'Hybrid BF', 'Hybrid BF with RF minimization', 'Hybrid BF with RF exhaustive search', 'Proposed', 'Location', 'Best');
lg.FontSize = 8;

savefig(['./results/ee_rate_with_eta_and_zeta_',num2str(Nt),'_',num2str(Lt),'_',num2str(Nr),'.fig'])
saveas(gcf,['./results/ee_rate_with_eta_and_zeta_',num2str(Nt),'_',num2str(Lt),'_',num2str(Nr),'.eps'],'epsc')

