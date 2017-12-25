clc;
clear all;
close all;

% calculate the SNR with Tx limitation and the noise parameter Tx power is
% 20dBm, and the noise is 3dB to -174dBm/Hz, so the noise db is
% -171dBm/Hz,and the pathloss is 120dB, so the receiver get -100dBm power maximum,
% assume the bandwidth BW=1MHZ, the noise -171dbm+10*lg(10^6)dB=-111dBm, so
% the SNR=-100+111=11dB, we can find the suitable SNR for the BER(bit error
% rate)<10^-3

% parasitics
Fc=2.4e8;   % carrier frequency
SNR=1:11;   % SNR range
BperS=1e6;  % Mbps value
N=1e4; % Number of bits
x_inp= randint(1,N);  % binary signal 0 or 1 % message to be transmitted                               
Tb=1/BperS; % bit period (second) =1/BperS   


% ========input signal with upsampling=====================================
x_bit=[]; 
nb=4; % upsampling rate
for n=1:1:N   % 
    if x_inp(n)==1;  
       x_bitt=ones(1,nb);
    else x_inp(n)==0;
        x_bitt=zeros(1,nb);
    end
    x_bit=[x_bit x_bitt];
end

t1=Tb/nb:Tb/nb:nb*N*(Tb/nb); % time of the signal 
f1=figure(1);
set(f1,'color',[1 1 1]);
subplot(2,1,1);
plot(t1,x_bit,'lineWidth',2);grid on;
axis([ 0 Tb*N -0.5 1.5]);
ylabel('Amplitude');
xlabel(' Time(sec)');
title('Input signal');


% BPSK Modulation
fi1=0; % carrier phase for bit 1
fi2=pi; % carrier phase for bit 0
t2=Tb/nb:Tb/nb:Tb;          
t2L=length(t2);
x_mod=[];
for (i=1:1:N)
    if (x_inp(i)==1)
        x_mod0=cos(2*pi*Fc*t2+fi1);%modulation signal 1
    else
        x_mod0=cos(2*pi*Fc*t2+fi2);%modulation signal 0
    end
    x_mod=[x_mod x_mod0];
end
t3=Tb/nb:Tb/nb:Tb*N;
subplot(2,1,2);
plot(t3,x_mod);
xlabel('Time(sec)');
ylabel('Amplitude');
title('Signal of BPSK modulation');


% =========transimission signal==========================================
x=x_mod;
% delay for 2us here data bit rate is 1M the time is 1 us and with 10
% upsampling, the t is 0.1us, so 2us delay is 20 bit shift
x_delay=[zeros(1,20) x_mod(1:length(x_mod)-20)];
x_mul=(x_mod+x_delay)*sqrt(1/2);% delay ISI signal

BER=[];
for i=1:11
    % ==========AWGN nosie====================================================
    n0=awgn(x,i);

    y0=x_mul+n0;
    % figure(4)
    % plot(t1,y0)
    % title('receiver signal with ISI');
    % ylabel('Amplitude');
    % xlabel('time(sec)');
    yFilt = conv(y0,ones(1,nb)); % convolution
    ySamp = yFilt(nb:nb:N*nb);  % sampling at time T=nb
    tnb=Tb:Tb:N*Tb;
    tfilt=Tb/nb:Tb/nb:length(yFilt)*Tb/nb;
    % figure(5)
    % subplot(2,1,1)
    % plot(tnb,ySamp)
    % title('downsampling convolution samples');
    % ylabel('Amplitude');
    % xlabel('time(sec)');
    % hold on
    % subplot(2,1,2)
    % plot(tfilt,yFilt)
    % title('receiver signal convolution for ISI demodulation');
    % ylabel('Amplitude');
    % xlabel('time(sec)');
    % making decision
    y_desion=[];
    for i=1:1:length(ySamp)
        if(ySamp(i)>0)
            y_d=1;
        else
            y_d=0;
        end
        y_desion=[y_desion y_d];
    end
    % BER
    error=0;
    for i=1:1:N
        if y_desion(i)~=x_inp(i)
            error=error+1;
        end
    end
    errorrate=error/N;
    BER=[BER errorrate];
end
figure(8)
plot(SNR,BER);
title('SNR-BER relation with upsampling rate at 4')
xlabel('SNR')
ylabel('BER')
% figure(20)
% plot(t1,n0);
% title('noise figure');
% ylabel('Amplitude');
% xlabel('time(sec)');


% % ==========receiver signal===============================================
% y=x+n0;
% figure(4)
% plot(t1,y)
% title('receiver signal');
% ylabel('Amplitude');
% xlabel('time(sec)');
% % =============demodulation================================================
% 
% % bypass filter for noise signal,give 1M BW
% [b,a] = butter(4,[2.4e8-1e6,2.4e8+1e6]/3e8);
% dataOut = filter(b,a,y);
% figure(5)
% plot(dataOut);
% title('after bypass filter signal');
% ylabel('Amplitude');
% xlabel('time(sec)');
% 
% % multi the carrier signal for next step, low pass filter
% y_1=[];
% for n=t2L:t2L:length(y)
%   t=Tb/nb:Tb/nb:Tb;
%   c=cos(2*pi*Fc*t); % carrier signal 
%   y_dem0=c.*y((n-(t2L-1)):n);
%   y_1=[y_1 y_dem0];
% end
% 
% % low pass filter
% Wc=2*1e8/Fc;                                          
% [b,a]=butter(4,Wc);
% y_dem = filter(b,a,y_1);
% figure(7)
% plot(y_dem);
% title('after low pass filter signal');
% ylabel('Amplitude');
% xlabel('time(sec)');
% 
% % making decision
% y_desion=[];
% for i=nb/2:nb:length(y_dem)
%     if(y_dem(i)>0)
%         y_d=1;
%     else
%         y_d=0;
%     end
%     y_desion=[y_desion y_d];
% end
% 
% x_out=y_desion; % output signal;
% 
% % BER
% error=0;
% for i=1:1:N
%     if x_out(i)~=x_inp(i)
%         error=error+1;
%     end
% end
% errorrate=error/N;


% 
% % ==============BER analysis==============================================
% % BER trace with the upsampling rate at 10, try different SNR to get BER
% BER=[];
% for i=1:11
%     n0=awgn(x,i);
%     
%     y=x+n0;
% 
% 
%     [b,a] = butter(4,[2.4e8-0.5e6,2.4e8+0.5e6]/3e8);
%     dataOut = filter(b,a,y);
% 
%     y_1=[];
%     for n=t2L:t2L:length(y)
%       t=Tb/nb:Tb/nb:Tb;
%       c=cos(2*pi*Fc*t); % carrier signal 
%       y_dem0=c.*y((n-(t2L-1)):n);
%       y_1=[y_1 y_dem0];
%     end
% 
% 
%     Wc=2*1e8/Fc;                                          
%     [b,a]=butter(4,Wc);
%     y_dem = filter(b,a,y_1);
% 
% 
%     y_desion=[];
%     for i=nb/2:nb:length(y_dem)
%         if(y_dem(i)>0)
%             y_d=1;
%         else
%             y_d=0;
%         end
%         y_desion=[y_desion y_d];
%     end
% 
%     x_out=y_desion; % output signal;
% 
%     % BER
%     error=0;
%     for i=1:1:N
%         if x_out(i)~=x_inp(i)
%             error=error+1;
%         end
%     end
%     errorrate=error/N;
%     BER=[BER errorrate];
% end
% 
% figure(8)
% plot(SNR,BER);
% title('SNR-BER relation with upsampling rate at 10')
% xlabel('SNR')
% ylabel('BER')
% 
% % =============upsampling rate analysis, with SNR=6========================
% nb_d=[2,4,6,8,10,20,40,100]; % different upsampling rate
% BER_nb=[];
% for i=1:length(nb_d)
%     nb=nb_d(i);
%     for n=1:1:N   
%         if x_inp(n)==1;  
%            x_bitt=ones(1,nb);
%         else x_inp(n)==0;
%             x_bitt=zeros(1,nb);
%         end
%         x_bit=[x_bit x_bitt];
%     end
% 
%     t1=Tb/nb:Tb/nb:nb*N*(Tb/nb); % time of the signal 
% 
% 
%     t2=Tb/nb:Tb/nb:Tb;          
%     t2L=length(t2);
%     x_mod=[];
%     for (i=1:1:N)
%         if (x_inp(i)==1)
%             x_mod0=cos(2*pi*Fc*t2+fi1);%modulation signal 1
%         else
%             x_mod0=cos(2*pi*Fc*t2+fi2);%modulation signal 0
%         end
%         x_mod=[x_mod x_mod0];
%     end
%     t3=Tb/nb:Tb/nb:Tb*N;
%     x=x_mod;
% 
%     n0=awgn(x,6);
% 
%     y=x+n0;
% 
% 
%     % bypass filter for noise signal
%     [b,a] = butter(4,[2.4e8-0.5e6,2.4e8+0.5e6]/3e8);
%     dataOut = filter(b,a,y);
% 
%     % multi the carrier signal for next step, low pass filter
%     y_1=[];
%     for n=t2L:t2L:length(y)
%       t=Tb/nb:Tb/nb:Tb;
%       c=cos(2*pi*Fc*t); % carrier signal 
%       y_dem0=c.*y((n-(t2L-1)):n);
%       y_1=[y_1 y_dem0];
%     end
% 
%     % low pass filter
%     Wc=2*1e8/Fc;                                          
%     [b,a]=butter(4,Wc);
%     y_dem = filter(b,a,y_1);
% 
%     % making decision
%     y_desion=[];
%     for i=nb/2:nb:length(y_dem)
%         if(y_dem(i)>0)
%             y_d=1;
%         else
%             y_d=0;
%         end
%         y_desion=[y_desion y_d];
%     end
% 
%     x_out=y_desion; % output signal;
% 
%     % BER
%     error=0;
%     for i=1:1:N
%         if x_out(i)~=x_inp(i)
%             error=error+1;
%         end
%     end
%     errorrate=error/N;
%     BER_nb=[BER_nb errorrate];
% end
% 
% figure(9)
% plot(nb_d,BER_nb);
% title('upsampling-BER relation with SNR=6')
% xlabel('upsampling rate')
% ylabel('BER')
