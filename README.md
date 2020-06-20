# SFND Radar Target Generation and Detection

<img src="https://github.com/epoc88/SensorFusion_RadarTargetGenerationAndDetection/blob/master/media/Radar_Project_Layout.png" width="700" height="400" />

---
#### 1. FMCW Waveform Design
Using the given system requirements, design a FMCW waveform. Find its Bandwidth (B), chirp time (Tchirp) and slope of the chirp.

```Matlab
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

CarrierFrequency= 77e9;             %
SweepTimeFactor = 5.5;

BandWith      = SpeedOfLight / (2 * RangeResolution);   % Bandwidth of the FMCW, Bsweep 
Tchirp = (SweepTimeFactor*2*MaxRange)/SpeedOfLight;     % Chirp Time of the FMCW
Slope  = BandWith/Tchirp;                               % Slope of the FMCW
```

#### 2. Simulation Loop
Simulate Target movement and calculate the beat or mixed signal for every timestamp.

```Matlab
TargetRange = 100;    
TargetVelocity = 50;  


Nd=128;          % number of chirps
Nr=1024;         % for length of time OR # of range cells

t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx = zeros(1,length(t)); %transmitted signal
Rx = zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t = zeros(1,length(t));
td = zeros(1,length(t));

for i=1:length(t)         
    
    r_t(i) = TargetRange + (TargetVelocity*t(i)); % range_covered
    td(i) = (2*r_t(i)) / SpeedOfLight; % time delay
    
    Tx(i) = cos(2 * pi * (CarrierFrequency * t(i) + (Slope * t(i) ^ 2 / 2)));
    Rx (i) = cos(2 * pi * (CarrierFrequency * (t(i) - td(i)) + (Slope * (t(i) - td(i)) ^ 2 / 2)));
    
    Mix(i) = Tx(i).*Rx(i);
    
end
```

#### 3. Range FFT (1st FFT)

Implement the Range FFT on the Beat or Mixed Signal and plot the result.

```Matlab
Mix = reshape(Mix, [Nr, Nd]);
SigFFT1 = fft(Mix, Nr) ./ Nr;
SigFFT1 = abs(SigFFT1); 

% Output of FFT is double sided signal, only one side of the spectrum is
% needed, thus we throw out half of the samples.
SigFFT1 = SigFFT1(1:Nr / 2);

%plotting the range
figure ('Name','Range from First FFT')
plot(SigFFT1);
axis ([0 200 0 1]);
xlabel('measured range');
```
<img src="https://github.com/epoc88/SensorFusion_RadarTargetGenerationAndDetection/blob/master/media/Range_1stFFT.jpg" width="700" height="400" />

#### 4. doppler FFT (2st FFT)

```Matlab
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
```
<img src="https://github.com/epoc88/SensorFusion_RadarTargetGenerationAndDetection/blob/master/media/RangeDopplerMap.jpg" width="700" height="400" />

#### 5. 2D CFAR
Implement the 2D CFAR process on the output of 2D FFT operation, i.e the Range Doppler Map.

Determine the number of Training cells for each dimension. Similarly, pick the number of guard cells.

```Matlab
Tr = 10;
Td = 8;
Gr = 4;
Gd = 4;
      
offset = 1.4;  % offset the threshold by SNR value in dB
```

Create a vector to store noise_level for each iteration on training cells，and  get the training Cells Num.

```Matlab
%Create a vector to store noise_level for each iteration on training cells
%noise_level = zeros(1,1);
noiseLevel = zeros(Nr/2-2*(Td+Gd),Nd-2*(Tr+Gr));
gridSize = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1);
trainingCellsNum = gridSize-(2*Gr+1)*(2*Gd+1);
CFARSig = zeros(size(RDM));
```

Slide the cell under test across the complete matrix. Make sure the CUT has margin for training and guard cells from the edges.
At every iteration, converting the value from logarithmic to linear using db2pow function, then sum the signal level within all the training cells, and obtain the mean value. then convert it back to db using pow2db. Adding the offset to it to determine the threshold. Then, comparing the signal under cut against this threshold. If the cut level > threshold assign it a value of 1, else equate it to 0.  Training, guard cells and offset are selected by increasing and decreasing to match the image shared in walkthrough.


```Matlab
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

```



<img src="https://github.com/epoc88/SensorFusion_RadarTargetGenerationAndDetection/blob/master/media/2DCFAR.jpg" width="700" height="400" />
