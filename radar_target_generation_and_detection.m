
clear; close all; clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

MaxRange = 200;
RangeResolution = 1;
MaxVelocity = 100;
SpeedOfLight = 3e8;

%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant

TargetRange = 100;    
TargetVelocity = 50;  
 


%% FMCW Waveform Generation

% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.


% Radar's carrier frequency
CarrierFrequency= 77e9;             %
SweepTimeFactor = 5.5;

BandWith      = SpeedOfLight / (2 * RangeResolution);   % Bandwidth of the FMCW, Bsweep 
Tchirp = (SweepTimeFactor*2*MaxRange)/SpeedOfLight;     % Chirp Time of the FMCW
Slope  = BandWith/Tchirp;                               % Slope of the FMCW
                                                          
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd = 128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr = 1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx = zeros(1,length(t)); %transmitted signal
Rx = zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t = zeros(1,length(t));
td = zeros(1,length(t));


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i=1:length(t)         
    
    r_t(i) = TargetRange + (TargetVelocity*t(i)); % range_covered
    td(i) = (2*r_t(i)) / SpeedOfLight; % time delay
    
    Tx(i) = cos(2 * pi * (CarrierFrequency * t(i) + (Slope * t(i) ^ 2 / 2)));
    Rx (i) = cos(2 * pi * (CarrierFrequency * (t(i) - td(i)) + (Slope * (t(i) - td(i)) ^ 2 / 2)));
    
    Mix(i) = Tx(i).*Rx(i);
    
end

%% RANGE MEASUREMENT


 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix, [Nr, Nd]);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
SigFFT1 = fft(Mix, Nr) ./ Nr;

% Take the absolute value of FFT output
SigFFT1 = abs(SigFFT1);  

 % *%TODO* :
% Output of FFT is double sided signal, only one side of the spectrum is
% needed, thus we throw out half of the samples.
SigFFT1 = SigFFT1(1:Nr / 2);


%plotting the range
figure ('Name','Range from First FFT')


 % plot FFT output 
plot(SigFFT1);
 
axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix = reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
SigFFT2 = fft2(Mix,Nr,Nd);

% Taking one side of signal from Range dimension.
SigFFT2 = SigFFT2(1:Nr/2,1:Nd);
SigFFT2 = fftshift (SigFFT2);
RDM = abs(SigFFT2);
RDM = 10*log10(RDM) ;

% leverage the surf function to plot the output of 2DFFT and to show axis in both dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr = 10;
Td = 8;

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr = 4;
Gd = 4;

% *%TODO* :
% offset the threshold by SNR value in dB
offset = 1.4;

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
%noise_level = zeros(1,1);
noiseLevel = zeros(Nr/2-2*(Td+Gd),Nd-2*(Tr+Gr));
gridSize = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1);
trainingCellsNum = gridSize-(2*Gr+1)*(2*Gd+1);

% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.
CFARSig = zeros(size(RDM));

   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR

for j=1:Nd-2*(Tr+Gr)
    for i=1:Nr/2-2*(Td+Gd)
        
        trainingCellsPatch = db2pow(RDM(i:i+2*(Td+Gd),j:j+2*(Gr+Tr)));
        trainingCellsPatch(Td+1:end-Td,Tr+1:end-Tr) = 0;
        
        noiseLevel(i,j) = pow2db(sum(sum(trainingCellsPatch))/trainingCellsNum);
        threshold = noiseLevel(i,j)*offset;
        
        if RDM(i+(Td+Gd),j+(Td+Gr))>threshold
            CFARSig(i+(Td+Gd),j+(Td+Gr)) = 1;
        else
            CFARSig(i+(Td+Gd),j+(Td+Gr)) = 0;
        end
           
    end
end



% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
 
RDM(RDM~=0 & RDM~=1) = 0;


% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
%figure,surf(doppler_axis,range_axis,'replace this with output');
figure('Name','CA-CFAR Filtered RDM'),surf(doppler_axis,range_axis,CFARSig);
colorbar;


 
 