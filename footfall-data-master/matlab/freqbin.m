function result=FREQBIN(name1,sr,slice_length,nchan,channo,switch,overlap)

% FREQBIN  result=FREQBIN(name,sr,slice_length,nchan,channo,switch,overlap)
%          The result vector contains rms 1/3 octave [4 Hz to half sample freq]
%          analysis of the vibration signal contained in the binary data file
%          NAME which contains data in IEE32 format. The rmq, PPV and overall
%          rms vibration levels are returned too.
%
%          result=[PPV rmq sig_len sr sl 1/3-Octave-RMS OverallRMS]
%
%          sig_len = length selected for analysing RMS data [s]
%          sr = sample rate [Hz]
%          sl = length of whole record [Hz]
%          slice_length = 2^n analysis slice (256, 512, 1024, 2048 etc)
%          nchan = no. of channel in file (e.g. 16)
%          channo = channel no. to be anlysed (base 0)
%          'switch' defines the way in which the programme operates.
%          		switch = 'a' is the automatic mode which analyses the event 
%          			between the 10dB down points [not working]
%          		switch =  'i' runs the programme interactively so that the 
%          			operator defines which part of the signal is analysed.
%          		switch=[t1 t2] will apply the analysis between times t1 and t2.
%
%          overlap should be set to 'y'.
%          The results are also written away to a .mat file "resbin.mat"
%          that is updated by each new analysis (i.e. if 25 files with 10 
%          channels are analyses, res_bin will have 250 lines of results.
%
%          N.B. jpeg versions of two key plots are stored to disc.
%
%
%          See also WRESULTS.M.
%          This file writes the file RES_BIN.MAT to the hard drive as a 
%          text format list of files and channels analysed and a text version of
%          the results.  This enables the data to be imported to spreadsheets.


%Programme Written By : Greer R,J 
%Updated : Subplot(211) section altered to include set scaleing of the plot
%          so that the data is not cutailed.
%          26/3/96, Conversion to Windows, IEE32 Input File, Multiple Channels
%	        Revised by R Greer 5 July 1999 to analyse data recorded by SNCF
%	        on Vouvray viaduct for CTRL Project
%          Revised 6 Nob 2000 to analyse data from new Arup CED acquisition system

% ******    Load data    ******
name2=[name1,'.bin'];
fid=fopen(name2);
DATA=fread(fid,inf,'float32');
fclose('all');


% ******    Initial Setup     ******
format long

figure(1);clf;
figure(2);clf;

disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp('       ******************************************')
disp('                                                 ')
disp('       Low Frequency Noise and Vibration Analysis')
disp('                                                 ')
disp('       Vertical Vibration Module: FREQBIN')
disp('                                                 ')
disp('       Written By : R,J Greer')
disp('                                                 ')
disp('       Version: 1.0')
disp('       ******************************************')


% ******    Initial parameters     ******
Header=DATA(1:5);% Remove header added by CED routine
DATA=DATA(6:length(DATA));
LTOT=length(DATA); % TOTAL NUMBER OF SAMPLES (ALL CHANNELS)
chans=length(channo);
ll=length(DATA)/nchan;   % TOTAL NUMBER OF SAMPLES (ONE CHANNEL)

if LTOT/nchan<=slice_length
   error('Data record is too short to use.  Choose smaller slice_length')
end

% ******   Reformat data into channels  *****
DATA=reshape(DATA,nchan,ll);


% ******   Extract first data channel   *****
data=DATA((channo(1)+1),:)';
name=[name1,'c',num2str(channo(1))];

% Remove D.C. Component
data=data-ones(ll,1)*mean(data);

time=(0:(ll-1))/sr;% TIME BASE VECTOR
t=1:16:ll;         % ONE IN 16 PLOT VECTOR
sl=max(time);     % TOTAL SAMPLE LENGTH
dmax=max(data);dmin=min(data);

figure(1)		
% ******    Define Evaluation of 10 dB down points     ******
if switch=='i'
   subplot(211);axis([0,max(time),1.2*dmin,1.2*dmax]);
     plot(time(t),data(t),'r')
	  xlabel('Data samples')
	  	ylabel('Un Cal Amplitude')
      titlestr=['Vibration Analysis Module for Event : ',name];
       title(titlestr)
		 hold on
       subplot(212);text(0,.9,'Please Use Mouse to Pick Start of Analysis Period');axis('off')
		 subplot(211);[x(1),crap]=ginput(1);plot([x(1) x(1)],[dmin dmax],'g--');
		 subplot(211),text(x(1),1.1*dmin,'Start')
       subplot(212);text(0,.7,'Please Use Mouse to Pick End of Analysis Period')
     	 subplot(211);[x(2),crap]=ginput(1);plot([x(2) x(2)],[dmin dmax],'g--');
		 subplot(211),text(x(2),1.1*dmin,'Finish')
       x=sr*x;
	    subplot(212);text(0,.5,'*** PRESS RETURN TO CONTINUE ***');
	  pause

elseif switch=='a'   
     % Evaluate 10dB down points
     x=find(data>=max(data)/3.16);
else
     x(1)=max(find(time<=switch(1)));
     x(2)=max(find(time<=switch(2)));
end

x1=min(x);x2=max(x);         % DEFINE TIMES OF 10 dB DOWN POINTS

sig_len=time(x2)-time(x1);  % SIGNAL LENGTH 
ld=length(data(x1:x2));			    % TOTAL NUMBER OF FREQ ANALYSIS POINTS

if ld<=slice_length*1.5
   error('Analysis period is too short to use.  Choose longer analysis period')
end

if ld>=LTOT/nchan-slice_length
   error('Analysis period is too long to use.  Choose a shorter analysis period')
end


% *****  Loop to analyse all channels  ******

for kk=1:chans
  name=[name1,'c',num2str(channo(kk))];
  data=DATA((channo(kk)+1),:)';

 % Remove D.C. Component
 data=data-ones(ll,1)*mean(data);

 dmax=max(data);dmin=min(data);


 figure(2);clf;
 % ******   PLOT SCALED TIME HISTORY      ******  
 subplot(211),axis([0,max(time),1.5*dmin,1.5*dmax]);% SET SCALING
 plot(time(t),data(t),'k')
 hold on
 subplot(211),plot([time(x1) time(x1)],1.6*[dmin dmax],'r');
 text(time(x1),1.2*dmax,'Start')
 subplot(211),plot([time(x2) time(x2)],1.6*[dmin dmax],'r');
 text(time(x2),1.2*dmax,'Finish')
 xlabel('Time [s]');ylabel('Velocity [mm/s]');
 titlestr=['Vibration Analysis Module for Event : ',name];
 title(titlestr)
 
 if overlap=='y'
    % ******     REARRANGE DATA IN TO SLICES WITH OVERLAP     ******
   N=floor(ld/slice_length)*2;
     if ld~=(N*(slice_length/2))
          vert1=data(x1:x1+(N*(slice_length/2))+((slice_length/2)-1));            % Truncate to slice_length block array.
     end
     Index=1:2*N*(slice_length/2);
     Index=reshape(Index,slice_length,N);
     Index=Index-ones(slice_length,1)*[0 (slice_length/2):(slice_length/2):(N-1)*(slice_length/2)];

  
    % *******   PLOT SLICE INFO TO TIME HISTORY     ******
    for iii=1:N
  	   	 yyy=0.9 + iii*(1.5-0.9)/N;
   	    plot([time(x1+Index(1,iii)) time(x1+Index(slice_length,iii))],yyy*[dmin dmin],'g')
	 		 plot([time(x1+Index(1,iii)) time(x1+Index(1,iii))],[(yyy-.1)*dmin yyy*dmin],'g')
	 		 plot([time(x1+Index(slice_length,iii)) time(x1+Index(slice_length,iii))],[0 yyy*dmin],'g:')
	%       text(time(x1+Index(1536,iii)),yyy*dmin,num2str(iii))
     end


	vert1=vert1(Index(:));
   vert1=reshape(vert1,slice_length,N);
 else
   % ******     REARRANGE DATA IN TO SLICES (NO OVERLAP)    *****
    vert1=data(x1:x2);
    l=length(vert1);
    N=floor(l/slice_length);% NUMBER OF SLICES
    if l~=(N*slice_length)
        vert1((N*slice_length)+1:l)=[];% TRUNCATE TO N*slice_length BLOCK VECTOR.
    end												  
    vert1=reshape(vert1,slice_length,N);% RESHAPE IN TO MATRIX/ARRAY.
   % *******   PLOT SLICE INFO TO TIME HISTORY     ******
   for iii=1:N
	   zzz=x1+iii*slice_length;
	   yyy=1.1 + iii*(1.5-1.1)/N;
      plot([time(zzz-slice_length) time(zzz)],yyy*[dmin dmin],'g')
	   plot([time(zzz-slice_length) time(zzz-slice_length)],[(yyy-.1)*dmin yyy*dmin],'g')
	   plot([time(zzz) time(zzz)],[0 yyy*dmin],'g:')
	   text(time(zzz-(slice_length/2)),yyy*dmin,num2str(iii))
   end
 end
 hold off

 % ************************************************************ 
 % ******     Beginning of Frequency Analysis Section     *****
 % ************************************************************

 % ******     Create hanning window and do FFTs.     ******
 w=hanning(slice_length)*ones(1,N);
 fdata=abs(fft(w.*vert1)/(slice_length/2));
    
 % ******     Generate frequency scale and prefered 1/3 octave band cutt offs
 freq=(0:slice_length/2)*(sr/slice_length);
 bands=[4 5 6.3 8 10 12.5 16 20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500]/(2^(1/6));
 f=find(freq<=1.1*max(bands));
 freq=freq(f);
 lf=length(f);
 fdata=fdata(f,:);
 fdata=fdata.^2;

max(10*log10(mean(fdata')/1e-12))


%figure(3);clf;mesh(flipud(10*log10(fdata)));
figure(3);clf;plot(freq,mean(10*log10(fdata'/1e-12)));
xlabel('Frequency [Hz]')
ylabel('Un Cal Mean RMS [dB]') 
 titlestr=['Analysis of Event : ',name];
 title(titlestr)
 
 eval(['print -djpeg90 ',name,'c']) %added by DMH to save spectrum plot (22/11/00)


 % ******    Evaluate Total and Percentage Energy in Each Band
 % ******    and Each Slice     ******************************
 
 for a=1:length(bands)
      index(a)=max(find(freq<=bands(a)));
 end
 totalenergy=sum(fdata((2:lf),:));% Total energy in each slice
 lb=length(index)-1;
 bndengy=zeros(lb,N);% Set up band energy array
 for a=1:lb
        [q,TT]=size(fdata((index(a)+1):index(a+1),:));
          if q==1,
              bndengy(a,:)=fdata(index(a)+1:index(a+1),:);
            elseif q==0,
              bndengy(a,:)=zeros(1,N);
            else
              bndengy(a,:)=sum(fdata(index(a)+1:index(a+1),:));
          end
 end
 
 % ****************************************************** 
 % ******     END of Frequency Analysis Section     *****
 % ******************************************************
  
 % ******     Calculate the RMQ, vrmsall, PPV and mean square values     ******
 rmq=(mean(vert1(:).^4))^0.25;
 VDVv=50.3/1000*rmq*sig_len^0.25;
 ms=mean(vert1.^2);
 vrmsall=sqrt(sum(vert1(:).^2)/(x2-x1));
 PPV=max([max(vert1(:)) -1*min(vert1(:))]); %edited by DMH to use cursor positons for rmq, vrmsall and ppv 22/11/00
 nl=0;% Digital Resolution 
 
 % ******     Evaluate band rms levels      ******
 ms=ones(lb,1)*ms;
 totalenergy=ones(lb,1)*totalenergy;
 ms=(ms.*(bndengy./totalenergy));
 bndrms=sqrt(mean(ms'));
waterfall=sqrt(ms);


 
 %******     Plot min max and mean 1/3 octave levels     ******
 hold off
 clear x y
 maxlevel=max(waterfall');
 minlevel=min(waterfall');
figure(2); subplot(212)
 
 urgh=find(maxlevel==0);maxlevel(urgh)=1e-5*ones(size(urgh));
 urgh=find(minlevel==0);minlevel(urgh)=1e-5*ones(size(urgh));
 urgh=find(bndrms==0);bndrms(urgh)=1e-5*ones(size(urgh));
 
 pdata1=20*log10(maxlevel(1:22)/1e-6);
 pdata2=20*log10(bndrms(1:22)/1e-6);
 %pdata3=20*log10(minlevel(1:22)/1e-6);

%size(pdata1)
%size(bands)
%pause
  semilogx(bands(1:22),pdata1,'r');axis([1,1000,max(pdata1)-60,max(pdata1)+15]);hold on
xlabel('Frequency [Hz]');
ylabel('RMS Band Level [dB re 1e-6 mm/s]')
grid

eval(['print -djpeg90 ',name,'a'])

figure(4);clf;
 [o1,o2]=bar(pdata1);fill(o1,o2,'r');axis([0,22,max(pdata1)-60,max(pdata1)+15]);hold on
 [o1,o2]=bar(pdata2);fill(o1,o2,'y');
 %[o1,o2]=bar(pdata3);fill(o1,o2,'b');
 titlestr=['Analysis of Event : ',name];
 title(titlestr)
 xlabel('1/3 O.B: 5=10Hz, 10=31.5Hz, 15=100Hz')
 ylabel('rms Vel Level [dB]')
 grid			   
 
 eval(['print -djpeg90 ',name,'b'])

 text(1,max(pdata1)+7,'KEY: Max-hold and Avg Spectra shown')
 
 % ******     SAVE DATA TO RESULTS FILE
 lamax=0;% Dummy for noise data
 if exist('res_bin.mat')~=2
    results=zeros(1,29);
    names='abcdefghhhh';
    save res_bin results names
 end
 
 load res_bin
 %bndrms=20*log10(bndrms/1e-6);
 %vrmsall=20*log10(vrmsall/1e-6);
 [n,m]=size(results);
 result=[PPV VDVv rmq sig_len sr sl bndrms(1:22) vrmsall] ;
 results(n+1,:)=result;
 if length(name)~=11
    for z=length(name):10
        name=[name,' '];
    end
 end
 names(n+1,:)=name;
 save res_bin results names 
 

 if switch=='i'
%  pause
 end

end;  % *** end kk loop ***
