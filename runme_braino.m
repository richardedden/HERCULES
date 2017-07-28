file='HB_9_2_raw_act.SDAT';
zerofillto=2048;
ppm=linspace(0,2000/128,zerofillto)-3.1325;
M.data=SDATreadMEGA(file,2048,128);
M.data_phase=M.data;
M.data_phase(1,1:2:end)=-M.data(1,2:2:end);

M.data=-M.data.*conj(repmat(M.data_phase(1,:),[2048 1]))./abs(repmat(M.data_phase(1,:),[2048 1]));
%remove T1 error from IWR
M.data(:,1:16:end)=M.data(:,1:16:end)*0.8;
figure(23)
plot(1:4096,real(fftshift(fft(M.data(1:1024,1:4:end),4096,1),1)),'b');
set(gca,'XLim',[1000 2000]);

M.data_avg=squeeze(sum(reshape(M.data,[2048 4 32]),3)),
 for ii=1:4
% %M.data(:,ii) = WaterFilter(M.data(:,ii),100,2000);
% 
 M.data_avg(:,ii) =waterremovalSVD(M.data_avg(:,ii), 2, 8, -0.08, 0.08, 0, 2048);
 end


%LB
LB=6;
window=repmat(exp(-(1:2048)/2000*pi*LB).',[1 4]);
%limit time domain
LIMIT_td=2048;
M.spec.ABCD=fftshift(fft(M.data_avg(1:LIMIT_td,:).*window(1:LIMIT_td,:),zerofillto,1));
%M.spec.ABCD=squeeze(sum(reshape(M.spec.raw,[2048 4 40]),3));
%Ph1
ph0=-24;
ph1=-100;
phase=repmat(exp(1i*pi*(ph1/180*(linspace(-0.5,0.5,zerofillto).')+ph0)),[1 4]);
M.spec.ABCD=M.spec.ABCD.*phase;
n_shift1=-0;
p_shift1=-3;
M.spec.ABCD(:,1)=circshift(M.spec.ABCD(:,1),[n_shift1 0])*exp(1i*pi/180*p_shift1);
n_shift2=0;
p_shift2=-3;
M.spec.ABCD(:,2)=circshift(M.spec.ABCD(:,2),[n_shift2 0])*exp(1i*pi/180*p_shift2);
n_shift3=0;
M.spec.ABCD(:,3)=circshift(M.spec.ABCD(:,3),[n_shift3 0]);
M.spec.H4=permute(hadamard(4)*permute(M.spec.ABCD,[2 1]),[2 1]);
%M.spec.H4(:,3)=-M.spec.H4(:,3);
M.spec.H4(:,1)=M.spec.H4(:,1);
M.spec.H4=circshift(M.spec.H4,[12 0]);
plot(ppm,M.spec.H4(:,:));
set(gca,'XLim',[0 4]);
%set(gca,'YLim',[-10 10]);
set(gca,'XDir','reverse');

fit_range = ppm >= 1.01 & ppm <= 4.5;
fit_range=circshift(fit_range,2);
ppm=ppm(fit_range);
sum(fit_range)

M.spec.H4=M.spec.H4(fit_range,:);


