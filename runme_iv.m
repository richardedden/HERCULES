file='hercmp_16_2_raw_act.SDAT';
ppm=linspace(0,2000/128,2048)-3.1325;
M.data=SDATreadMEGA(file,2048,160);


M.data=M.data.*conj(repmat(M.data(1,:),[2048 1]))./abs(repmat(M.data(1,:),[2048 1]));

M.data_avg=squeeze(sum(reshape(M.data,[2048 4 40]),3)),
% for ii=1:4
% %M.data(:,ii) = WaterFilter(M.data(:,ii),100,2000);
% 
% M.data_avg(:,ii) =waterremovalSVD(M.data_avg(:,ii), 2, 8, -0.08, 0.08, 0, 2048);
% end


%LB
LB=0;
window=repmat(exp(-(1:2048)/2000*pi*LB).',[1 4]);
%limit time domain
LIMIT_td=2048;
M.spec.ABCD=fftshift(fft(M.data_avg(1:LIMIT_td,:).*window(1:LIMIT_td,:),2048,1));
%M.spec.ABCD=squeeze(sum(reshape(M.spec.raw,[2048 4 40]),3));
%Ph1
ph1=-180;
phase=repmat(exp(1i*pi*ph1/360*(linspace(-0.5,0.5,2048).')),[1 4]);
M.spec.ABCD=M.spec.ABCD.*phase;
M.spec.H4=permute(hadamard(4)*permute(M.spec.ABCD,[2 1]),[2 1]);
%M.spec.H4(:,3)=-M.spec.H4(:,3);
M.spec.H4(:,1)=M.spec.H4(:,1);
plot(ppm,M.spec.H4);
set(gca,'XLim',[0 5]);
set(gca,'XDir','reverse');

fit_range = ppm >= 1.01 & ppm <= 4.5;
fit_range=circshift(fit_range,2);
ppm=ppm(fit_range);
sum(fit_range)

M.spec.H4=M.spec.H4(fit_range,:);


