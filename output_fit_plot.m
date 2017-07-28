%start off with the final plort of fit-IV_2.M
hold off
figure(2)

x_shift=500;
x_display=[(1:size(Sim_Spectrum_display,1)).' (1:size(Sim_Spectrum_display,1)).'+x_shift (1:size(Sim_Spectrum_display,1)).'+x_shift*2 (1:size(Sim_Spectrum_display,1)).'+x_shift*3 ];
size(x_display);
plot(x_display,Sim_Spectrum_display,'k',x_display,Fit_Spectrum_display,'r',x_display,spline_vault_display,'g')

y_shift=1E5;

hold on

Names{1}='Asc';
Names{2}='Asp';
Names{3}='Cr';
Names{4}='GABA';
Names{5}='GPC';
Names{6}='GSH';
Names{7}='Gln';
Names{8}=ANames{13}(3:end);
Names{9}=ANames{14}(3:end);
Names{10}=ANames{15}(3:end);
Names{11}=ANames{16}(3:end);
Names{12}='NAA';
Names{13}='NAAG';
Names{14}='PCho1';
Names{15}='PCho2';

Names{16}='PCr1';;
Names{17}='PCr2';;
Names{18}='Scyllo';;
Names{19}='2HG';


for jj=1:19
    Fitted_weights_dummy=Fitted_weights(1,:);
    Fitted_weights_dummy(20:27)=0;
    Fitted_weights_dummy(30)=0;
    Fitted_weights_dummy(jj)=Fitted_weights_dummy(jj)*5;
     if jj>1
         Fitted_weights_dummy(1:(jj-1))=0;
     end
     if jj<19
        Fitted_weights_dummy((jj+1):19)=0;
    end
    
plot(x_display,reshape(Hercules_model(Fitted_weights_dummy,BasisSet3),[457, 4])-y_shift*10-y_shift*jj*10,'k');
text(420,-y_shift*(jj+0.6)*10,Names{jj})
text(950,-y_shift*(jj+0.7)*10,num2str(Fitted_weights_dummy(jj)/6,2))
end
