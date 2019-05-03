
                           %%%%%%%%%%%%%%%%%%%%
                           %%% Stat_Builder %%%
                           %%%%%%%%%%%%%%%%%%%%
                           
% Written by Nijil Lal and Biveen Shajilal to build the statistics of
% Poissonian, super-Poissonian and sub-Poissonian light sources.

%% --------------------------- Reading Data --------------------------- %%

    % A usual oscilloscope records the waveworm with the timestamp in the
    % first column and the input channels in the onward columns.
    
data = csvread('/home/user/data.csv');

A = data(:,[1 2]);  % channel 1
B = data(:,[1 3]);  % channel 2

%% -------------------- Conversion to Binary Signals -------------------%% 

    % The recorded time series from oscilloscope is converted into a series
    % of binaries (1-event, 0-no event). To isolate the events from noise
    % signals and possible after pulses, a threshold (thresh) is determined
    % which is above the noise level and any undesired afterpulses. In case
    % of an event signal that crosses the threshold, 1 is assigned to the
    % binary series at the onset of the signal - defined by signal_pos,
    % whpse value is much smaller than the threshold.
    % Here, they are chosen as 0.4 and 0.1 respectively.

N = length (data);
    thresh = 0.4;
    signal_pos = 0.1;

A_binary(:,1) = A(:,1);
B_binary(:,1) = A(:,1);
AB_binary(:,1)= A(:,1);
AB_binary(1:N,2)=0;


for i=1:N
    
    if A(i,2)>thresh
        j=i;
         while A(j,2)>signal_pos
             A_binary(j,2)=0;
             j=j-1;
         end
         A_binary(j,2)=1;  
    else
        A_binary(i,2)=0;
    end

end


for i=1:N
    
    if B(i,2)>thresh
        j=i;
         while B(j,2)>signal_pos
             B_binary(j,2)=0;
             j=j-1;
         end
         B_binary(j,2)=1;  
    else
        B_binary(i,2)=0;
    end

end

% The coincidence time series is built from the binary series, A and B.

tc = 1; 

% tc x resolution (of the recorded time series) gives the width of the
% coincidence window.

for i=1:N
    
    if A_binary(i,2)==1
        for k=i:i+tc     
            if B_binary(k,2)==1
                AB_binary(k,2)=1;
            end
        end
    end
end

%% --------------- Determining the size of the time bin --------------- %%
   % Bin size is defined for rach series for a mean photon number = 1/bin.

nA = sum(A_binary(:,2));
nB = sum(B_binary(:,2));
nAB = sum(AB_binary(:,2));
% Total number of photons corresponding to each series.

binA = round(N/nA); binB = round(N/nB); binAB = round(N/nAB);
% Average photon number per bin is defined by the size of bin.

MaxA = N-binA; MaxB = N-binB; MaxAB = N-binAB;
% Upper limit while slicing the time series into bins.

%% --------------------- Building The Statistics --------------------- %%

% Arm A statistics
j=1;
for i=1 : binA : MaxA

A_count=0;
    for k=i:(i+binA)
        
        if A_binary(k,2)==1
        A_count=A_count+1;
        end
        
    end
    
    A_stat(j,1)=A_count;
    j=j+1;    
end

% Arm B statistics
j=1;
for i=1 : binB : MaxB
    
B_count=0;    
    for k=i:(i+binB)
        if B_binary(k,2)==1
        B_count=B_count+1;
        end
        
    end
    
    B_stat(j,1)=B_count;
    j=j+1;
end

% Heralded single photon statistics
j=1;
for i=1 : binAB : MaxAB

AB_count=0;
    for k=i:(i+binAB)
        if AB_binary(k,2)==1
        AB_count=AB_count+1;
        end
        
    end
    
    AB_stat(j,1)=AB_count;
    j=j+1;
end

%% ---------------------- Building the Histogram ---------------------- %%

% h is the number of photons,n, in each bin.

for h=1:20
    
    U(h,1)=sum(A_stat==h-1);
    V(h,1)=sum(B_stat==h-1);
    W(h,1)=sum(AB_stat==h-1);
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% ************************** ---- End ---- ************************** %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% %%