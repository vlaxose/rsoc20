clear
clc

users_range = [1 5 10 15 20];
Nt=32;Lt=8;Nr=8;K=5;

totalMCrealizations = 500;
rate = zeros(totalMCrealizations, length(users_range), 6);
power = zeros(totalMCrealizations, length(users_range), 6);
noise_variance = sqrt(10^(-5/10));

for user_index = 1:length(users_range)
    K = users_range(user_index);
    for r=1:totalMCrealizations
        [rate(r, user_index, :), power(r, user_index, :)] = systemModel(Nt, Nr, Lt, K, noise_variance);
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
p1=plot(users_range, mean_rate(:, 1)); hold on;
set(p1,'LineWidth',1, 'LineStyle', '-', 'Color', 'Black');
p2=plot(users_range, mean_rate(:, 2)); hold on;
set(p2,'LineWidth',1, 'LineStyle', '--', 'Color', 'Black');
p3=plot(users_range, mean_rate(:, 3)); hold on;
set(p3,'LineWidth',1, 'LineStyle', '-.', 'Color', 'Black');
p4=plot(users_range, mean_rate(:, 4)); hold on;
set(p4,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Green', 'Marker', 'o');
p5=plot(users_range, mean_rate(:, 5)); hold on;
set(p5,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Blue', 'Marker', 's');
p6=plot(users_range, mean_rate(:, 6)); hold on;
set(p6,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Red', 'Marker', 'h');
grid on;
xlabel('Number of Users', 'FontSize', 11)
ylabel({'Rate', '(Mbits/sec)'}, 'FontSize', 11)
lg = legend('Digital BF', 'Analog BF', 'Hybrid BF', 'Hybrid BF with RF minimization', 'Hybrid BF with RF exhaustive search', 'Proposed', 'Location', 'Best');
lg.FontSize = 8;

savefig(['./results/rate_users_',num2str(Nt),'_',num2str(Lt),'_',num2str(Nr),'.fig'])
saveas(gcf,['./results/rate_users_',num2str(Nt),'_',num2str(Lt),'_',num2str(Nr),'.eps'],'epsc')

figure;
mean_ee(:, 1) = mean_rate(:, 1)./mean_power(:, 1);
p1=plot(users_range, mean_ee(:, 1)); hold on;
set(p1,'LineWidth',1, 'LineStyle', '-', 'Color', 'Black');
mean_ee(:, 2) = mean_rate(:, 2)./mean_power(:, 2);
p2=plot(users_range, mean_ee(:, 2)); hold on;
set(p2,'LineWidth',1, 'LineStyle', '--', 'Color', 'Black');
mean_ee(:, 3) = mean_rate(:, 3)./mean_power(:, 3);
p3=plot(users_range, mean_ee(:, 3)); hold on;
set(p3,'LineWidth',1, 'LineStyle', '-.', 'Color', 'Black');
mean_ee(:, 4) = mean_rate(:, 4)./mean_power(:, 4);
p4=plot(users_range, mean_ee(:, 4)); hold on;
set(p4,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Green', 'Marker', 'o');
mean_ee(:, 5) = mean_rate(:, 5)./mean_power(:, 5);
p5=plot(users_range, mean_ee(:, 5)); hold on;
set(p5,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Blue', 'Marker', 's');
mean_ee(:, 6) = mean_rate(:, 6)./mean_power(:, 6);
p6=plot(users_range, mean_ee(:, 6)); hold on;
set(p6,'LineWidth',1.5, 'LineStyle', '-', 'Color', 'Red', 'Marker', 'h');
grid on;
xlabel('Number of Users', 'FontSize', 11)
ylabel({'Energy Efficiency', '(Mbits/Joule)'}, 'FontSize', 11)
lg = legend('Digital BF', 'Analog BF', 'Hybrid BF', 'Hybrid BF with RF minimization', 'Hybrid BF with RF exhaustive search', 'Proposed', 'Location', 'Best');
lg.FontSize = 8;

savefig(['./results/ee_users_',num2str(Nt),'_',num2str(Lt),'_',num2str(Nr),'.fig'])
saveas(gcf,['./results/ee_users_',num2str(Nt),'_',num2str(Lt),'_',num2str(Nr),'.eps'],'epsc')

