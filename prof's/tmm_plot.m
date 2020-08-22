 %% plot
 % display

% hold on
% titl=append('cavity width=',sprintf('%0.5g', cw),'nm')
% subplot(2,1,1), plot(lambda_a,Re_S), title(titl), legend('R')
% ylim([0,1])
% set(gca,'fontsize', 16);
% subplot(2,1,2), plot(lambda_a,Abs_S), legend('A')
% ylim([0,1])
% set(gca,'fontsize', 16);
% xlabel('wavelength')



% figure
% subplot(3,1,1), plot(lambda_a,Abs_S,lambda_a,Abs_P), title('Absorption Coeff'), legend('S','P')
% ylim([0,1])
% set(gca,'fontsize', 14);
% subplot(3,1,2), plot(lambda_a,Tr_S,lambda_a,Tr_P), title('Transmission Coeff'), legend('S','P')
% ylim([0,1])
% set(gca,'fontsize', 14);
% subplot(3,1,3), plot(lambda_a,Re_S,lambda_a,Re_P), title('Reflection Coeff'), legend('S','P')
% ylim([0,1])
% set(gca,'fontsize', 14);
% sgtitle('fun funky stuff')

% figure
% lambda_eV = (1240./(lambda_a/10^-9)).';
% eV_spaced = (linspace(lambda_eV(size(lambda_a,2)),lambda_eV(1),1000)).';
% Re_S_spaced = (interp1(lambda_eV,Re_S,eV_spaced));
% 
% imagesc(theta0degrees_a,eV_spaced, Re_S_spaced)
% axis xy
% xlabel('angle')
% ylabel('energy')
% ylim([1.5,3])
% set(gca,'fontsize', 16);

figure
%a0 and a'
plot(lambda_a,perovskiteinfo(1).A,lambda_a,perovskiteinfo(2).A), legend('A0','A1'), title('Absorption of Perovskite with 1% contraction'), xlabel('wavelength')
ylim([0 1])
figure
%r0 and r'
plot(lambda_a,finalinfo(1).R,lambda_a,finalinfo(3).R), legend('R0','R1'), title('reflection of cavity with 1% contraction'), xlabel('wavelength')
figure
% dA and DR
plot(lambda_a,perovskiteinfo(2).dA,lambda_a,finalinfo(3).dR), legend('dA','dR'), title('dA vs dR'), xlabel('wavelength')

% figure
% %for multiple dRs
% for c=3:1:3
%     hold on
%     plot(lambda_a,finalinfo(c).dR), title('dR')
%     %plot(lambda_a,perovskiteinfo(c).dA), title('dA')
%     %plot(lambda_a,finalinfo(c).R)
%     %plot(lambda_a,perovskiteinfo(c).A)
%     %imagesc(theta0degrees_a,lambda_a, finalinfo(c).dR)
% 
% end
