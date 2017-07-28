hold off
figure(3)
Fitted_weights(ii,12)=6;
Fitted_weights(ii,13)=0;


Fit_Spectrum2=reshape(Hercules_model(Fitted_weights(ii,:),BasisSet3),[length(ppm),4])+squeeze(sum(spline_vault,3));
Fit_Spectrum_display=Fit_Spectrum2;
Fit_Spectrum_display(:,2:4)=Fit_Spectrum_display(:,2:4)*5;

plot(x_display,Sim_Spectrum_display,'k',x_display,Fit_Spectrum_display,'r',x_display,spline_vault_display,'g')



% hold on
% 
% 
% 
% for jj=1:19
%     Fitted_weights_dummy=Fitted_weights(1,:);
%     Fitted_weights_dummy(20:27)=0;
%     Fitted_weights_dummy(30)=0;
%     Fitted_weights_dummy(jj)=Fitted_weights_dummy(jj)*5;
%      if jj>1
%          Fitted_weights_dummy(1:(jj-1))=0;
%      end
%      if jj<19
%         Fitted_weights_dummy((jj+1):19)=0;
%     end
%     
% plot(x_display,reshape(Hercules_model(Fitted_weights_dummy,BasisSet3),[457, 4])-y_shift*10-y_shift*jj*10,'k');
% text(420,-y_shift*(jj+0.6)*10,Names{jj})
% text(950,-y_shift*(jj+0.7)*10,num2str(Fitted_weights_dummy(jj)/3,2))
% end
