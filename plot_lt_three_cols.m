clear
clc

lt_range = [2 4 8 10 12];
snr=15;
noise_variance = sqrt(10^(-snr/10));    

Nt=32;Nr=8;K=8;Lmax=15;

totalMCrealizations = 1000;
rate = zeros(totalMCrealizations, length(lt_range), 6);
power = zeros(totalMCrealizations, length(lt_range), 6);

for lt_index = 1:length(lt_range)
    Lt = lt_range(lt_index);
    for r=1:totalMCrealizations
        [rate(r, lt_index, :), power(r, lt_index, :)] = systemModel(Nt, Nr, Lt, K, noise_variance,Lmax);
    end
end
mean_rate = squeeze(sum(rate, 1))/totalMCrealizations;
mean_power = squeeze(sum(power, 1))/totalMCrealizations;    

mean_ee = zeros(size(mean_rate));

figure;
subplot(1,3,1);
p1=plot(lt_range, mean_rate(:, 1)); hold on;
set(p1,'LineWidth',1, 'LineStyle', '-', 'Color', 'Black');
p2=plot(lt_range, mean_rate(:, 2)); hold on;
set(p2,'LineWidth',1, 'LineStyle', '--', 'Color', 'Black');
p3=plot(lt_range, mean_rate(:, 3)); hold on;
set(p3,'LineWidth',1, 'LineStyle', '-.', 'Color', 'Black');
p4=plot(lt_range, mean_rate(:, 4)); hold on;
set(p4,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Green', 'Marker', 'o');
p5=plot(lt_range, mean_rate(:, 5)); hold on;
set(p5,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Blue', 'Marker', 's');
p6=plot(lt_range, mean_rate(:, 6)); hold on;
set(p6,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Red', 'Marker', 'h');
grid on;
xlabel('Number of RF chains', 'FontSize', 11)
ylabel({'Spectral Efficiency', '(Mbits/sec)'}, 'FontSize', 11)

subplot(1,3,2);
p1=plot(lt_range, mean_power(:, 1)); hold on;
set(p1,'LineWidth',1, 'LineStyle', '-', 'Color', 'Black');
p2=plot(lt_range, mean_power(:, 2)); hold on;
set(p2,'LineWidth',1, 'LineStyle', '--', 'Color', 'Black');
p3=plot(lt_range, mean_power(:, 3)); hold on;
set(p3,'LineWidth',1, 'LineStyle', '-.', 'Color', 'Black');
p4=plot(lt_range, mean_power(:, 4)); hold on;
set(p4,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Green', 'Marker', 'o');
p5=plot(lt_range, mean_power(:, 5)); hold on;
set(p5,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Blue', 'Marker', 's');
p6=plot(lt_range, mean_power(:, 6)); hold on;
set(p6,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Red', 'Marker', 'h');
grid on;
xlabel('Number of RF chains', 'FontSize', 11)
ylabel({'Power', '(Joules/sec)'}, 'FontSize', 11)

subplot(1,3,3);

mean_ee(:, 1) = mean_rate(:, 1)./mean_power(:, 1);
p1=plot(lt_range, mean_ee(:, 1)); hold on;
set(p1,'LineWidth',1, 'LineStyle', '-', 'Color', 'Black');
mean_ee(:, 2) = mean_rate(:, 2)./mean_power(:, 2);
p2=plot(lt_range, mean_ee(:, 2)); hold on;
set(p2,'LineWidth',1, 'LineStyle', '--', 'Color', 'Black');
mean_ee(:, 3) = mean_rate(:, 3)./mean_power(:, 3);
p3=plot(lt_range, mean_ee(:, 3)); hold on;
set(p3,'LineWidth',1, 'LineStyle', '-.', 'Color', 'Black');
mean_ee(:, 4) = mean_rate(:, 4)./mean_power(:, 4);
p4=plot(lt_range, mean_ee(:, 4)); hold on;
set(p4,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Green', 'Marker', 'o');
mean_ee(:, 5) = mean_rate(:, 5)./mean_power(:, 5);
p5=plot(lt_range, mean_ee(:, 5)); hold on;
set(p5,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Blue', 'Marker', 's');
mean_ee(:, 6) = mean_rate(:, 6)./mean_power(:, 6);
p6=plot(lt_range, mean_ee(:, 6)); hold on;
set(p6,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Red', 'Marker', 'h');
grid on;
xlabel('Number of RF chains', 'FontSize', 11)
ylabel({'Energy Efficiency', '(Mbits/Joule)'}, 'FontSize', 11)

lg = legend('Digital BF', 'Analog BF', 'Hybrid BF', 'Hybrid BF with RF minimization', 'Hybrid BF with RF exhaustive search', 'Proposed');
lg.FontSize = 8;

savefig(['./results/ee_lt_three_cols_',num2str(Nt),'_',num2str(Lt),'_',num2str(Nr),'.fig'])
saveas(gcf,['./results/ee_lt_three_cols_',num2str(Nt),'_',num2str(Lt),'_',num2str(Nr),'.eps'],'epsc')
