% for replicating figures in xenon acceleration paper
close all;
clear

load("simresults.mat")


%% CFA
figure;tiledlayout(2,3);

plot_sim(io,B1_range,Noise_range)
plot_sim(io_pf,B1_range,Noise_range)

metric={"MSE"};
h1=figure('Name',sprintf('CFA',string(metric)));

compareMetrics(B1_range, Noise_range, metric,io,io_pf);%title('Inside-out')
legend('High noise','High noise (under-sampled)','Medium noise', 'Medium noise (under-sampled)','Low noise', 'Low noise (under-sampled)','Location','southoutside')
%% RFA
figure;tiledlayout(4,3);

plot_sim(io_ramped,B1_range,Noise_range)
plot_sim(io_ramped_pf,B1_range,Noise_range)

plot_sim(oi_ramped,B1_range,Noise_range)
plot_sim(oi_ramped_pf,B1_range,Noise_range)

metric={"MSE"};
%inout vs outin
h2=figure('Name',sprintf('IO vs OI (%s)',string(metric)));
tiledlayout(2,1)

compareMetrics(B1_range, Noise_range, metric,io_ramped,io_ramped_pf);title('Inside-out')
compareMetrics(B1_range, Noise_range, metric,oi_ramped,oi_ramped_pf);title('Outside-in')
legend('High noise','High noise (under-sampled)','Medium noise', 'Medium noise (under-sampled)','Low noise', 'Low noise (under-sampled)','Location','southoutside')

%% plot mean estimates
function plot_sim(struct,B1_range, Noise_range)
[X,Y] = meshgrid(B1_range, Noise_range);

nexttile;
surf(X, Y, mean(struct.MSE,3))
xlabel('\DeltaB_{1}^{+}');
ylabel('\sigma_{k}^2'); yticks(10.^[-2,0,2]);set(gca,'YScale','Log')
zlabel('MSE')
axis square;

nexttile;
surf(X, Y, mean(struct.PSNR,3))
xlabel('\DeltaB_{1}^{+}');
ylabel('\sigma_{k}^2'); yticks(10.^[-2,0,2]);set(gca,'YScale','Log')
zlabel('PSNR');set(gca, 'ZDir', 'reverse')
axis square;

nexttile;
surf(X, Y, mean(struct.SSIM,3))
xlabel('\DeltaB_{1}^{+}');
ylabel('\sigma_{k}^2'); yticks(10.^[-2,0,2]);set(gca,'YScale','Log')
zlabel('SSIM');set(gca, 'ZDir', 'reverse')
axis square;

end

%% compare
function compareMetrics(B1_range, Noise_range, metric,varargin)

colors=[213,94,0;0,0,0,;0,114,178]./256;
lines={'-','--'};

noise_ind=[24,20,6];
nexttile;hold on
for i=1:length(noise_ind)
    for j=1:numel(varargin)
    data=mean(varargin{j}.(string(metric)),3);
        plot(B1_range,data(noise_ind(i),:),'Linewidth',2,'LineStyle',lines{j},'Color',colors(i,:));
    end
end
axis tight;axis square
xlabel('\DeltaB_{1}^{+}');
ylabel('MSE');
ylim([0 0.06]);
% B1_ind=[6,16,26];
% nexttile;hold on
% for j=1:numel(varargin)
%     data=mean(varargin{j}.(string(metric)),3);
%     for i=1:length(B1_ind)
%         plot(Noise_range,data(:,B1_ind(i)),'Linewidth',2,'LineStyle',lines{j},'Color',colors(i,:));
%     end
% end
% axis square
% xlabel('Noise');
% % set(gca,'XScale','Log')
% ylabel('MSE');

end

function compareMetricsMap(B1_range, Noise_range, metric,case1,case2)
[cMap] = RGcMap(1);

data=mean(case1.(string(metric)),3)-mean(case2.(string(metric)),3);
nexttile;imagesc(B1_range, Noise_range, data);
clim([-1e-2 1e-2])
axis xy
axis square
xlabel('\DeltaB_{1}^{+}');
ylabel('Noise [dB]');
% set(gca,'YScale','Log')

colormap(cMap)
colorbar
end