% display
% 
% subplot(2,1,1), plot(lambda_a,Re_S,lambda_a,Re_P), title('Reflection Coeff'), legend('S','P')
% ylim([0,1])
% set(gca,'fontsize', 16);
% subplot(2,1,2), plot(lambda_a,Tr_S,lambda_a,Tr_P), title('Transmission Coeff'), legend('S','P')
% ylim([0,1])
% set(gca,'fontsize', 16);

subplot(3,1,1), plot(lambda_a,Re_S,lambda_a,Re_P), title('Reflection Coeff'), legend('S','P')
ylim([0,1])
set(gca,'fontsize', 14);
subplot(3,1,2), plot(lambda_a,Tr_S,lambda_a,Tr_P), title('Transmission Coeff'), legend('S','P')
ylim([0,1])
set(gca,'fontsize', 14);
subplot(3,1,3), plot(lambda_a,Abs_S,lambda_a,Abs_P), title('Absorption Coeff'), legend('S','P')
ylim([0,1])
set(gca,'fontsize', 14);

% lambda_eV = 1240./(lambda_a/10^-9);
% imagesc(theta0degrees_a,lambda_eV, Re_P)
% axis xy
% xlabel('angle')
% ylabel('energy')
% ylim([1.5,3.5])
% set(gca,'fontsize', 16);
