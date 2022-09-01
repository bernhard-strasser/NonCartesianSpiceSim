% Simulate Phantom Data



%% Create Phantom Masks


% Load water mask
load('./InputMeasAndLogData/SpiralSpiceSim_WaterMask2.mat')

% % For comparing the mask to the signal:
% load('../InputMeasAndLogData/B0Map.mat','B0Map_MagSoS')


SimPhantom.Masks.Vial{1} = mask;    % Water-Background Vial
clear mask;


% Create Vial Masks
% Vial 2
Vial = EllipticalFilter(ones(size(SimPhantom.Masks.Vial{1})),[1 2],[1 1 1 8.5],1);
Vial = circshift(Vial,-[65-53 65-51]);
SimPhantom.Masks.Vial{2} = Vial;

% Vial 3
Vial = EllipticalFilter(1*ones(size(SimPhantom.Masks.Vial{1})),[1 2],[1 1 1 3.5],1);
Vial = circshift(Vial,-[65-57 65-87]);
SimPhantom.Masks.Vial{3} = Vial;

% Vial 4
Vial = EllipticalFilter(1*ones(size(SimPhantom.Masks.Vial{1})),[1 2],[1 1 1 4.5],1);
Vial = circshift(Vial,-[65-43 65-67]);
SimPhantom.Masks.Vial{4} = Vial;

% Vial 5
Vial = EllipticalFilter(1*ones(size(SimPhantom.Masks.Vial{1})),[1 2],[1 1 1 3.5],1);
Vial = circshift(Vial,-[65-57 65-70]);
SimPhantom.Masks.Vial{5} = Vial;

% Vial 6
Vial = EllipticalFilter(1*ones(size(SimPhantom.Masks.Vial{1})),[1 2],[1 1 1 5.5],1);
Vial = circshift(Vial,-[65-77 65-76]);
SimPhantom.Masks.Vial{6} = Vial;

% Vial 7
Vial = EllipticalFilter(1*ones(size(SimPhantom.Masks.Vial{1})),[1 2],[1 1 1 8.5],1);
Vial = circshift(Vial,-[65-76 65-54]);
SimPhantom.Masks.Vial{7} = Vial;

clear Vial;


% imresize if necessary
for i = 1:numel(SimPhantom.Masks.Vial)
    SimPhantom.Masks.Vial{i} = imresize(SimPhantom.Masks.Vial{i},ParSim.MatSize(1:2),'nearest');
end



Masky = SimPhantom.Masks.Vial{2};
for i = 3:7
    Masky = Masky + (i-1)*SimPhantom.Masks.Vial{i};
end
SimPhantom.Masks.VialMask = logical(Masky);


if(~quiet_flag)
   Masky = Masky + SimPhantom.Masks.Vial{1};
   figure; imagesc(Masky); title('Phantom Masks')
   clear Masky
end


%% Create Spatial Metabolite Maps

SimPhantom.MetaboConc.Water = [24 10 10 10 10 10 10];          % No water for now
SimPhantom.MetaboConc.NAA = [0 1 2 3 4 5 6];
SimPhantom.MetaboConc.Cr = [0 6 5 4 3 2 1];
SimPhantom.MetaboConc.Cho = [0 3 1 4 6 2 5];
SimPhantom.MetaboConc.Glu = [0 6 1 3 5 2 4];       

% % Simulate water only
% SimPhantom.MetaboConc.Water = [24 1 2 3 4 5 6];          
% SimPhantom.MetaboConc.NAA = [0 0 0 0 0 0 0];
% SimPhantom.MetaboConc.Cr = [0 0 0 0 0 0 0];
% SimPhantom.MetaboConc.Cho = [0 0 0 0 0 0 0];
% SimPhantom.MetaboConc.Glu = [0 0 0 0 0 0 0];     



for CurField = transpose(fields(SimPhantom.MetaboConc))
   
    SimPhantom.MetaboMaps.(CurField{1}) = zeros(size(SimPhantom.Masks.Vial{1}));
    for VialNo = 1:numel(SimPhantom.Masks.Vial)
        SimPhantom.MetaboMaps.(CurField{1}) = SimPhantom.MetaboMaps.(CurField{1}) + ...
        SimPhantom.MetaboConc.(CurField{1})(VialNo) * SimPhantom.Masks.Vial{VialNo};
    end
    
end





%% Create Different Spectra for Different Vials & Water-Background

% Simulate FIDs

Dt = ParSim.SpecDwellTimes / 10^9;
% [FID,Spectrum]=Simulate_FID_Spectra(Chemshift,DeltaFrequency,phase0,AcqDelay,T2,S_0,SNR,dwelltime,vecSize,LarmorFreq,SmoothFIDSpan)
[SimPhantom.FIDs.Water, ppm] = Simulate_FID_Spectra(4.65,4.65,0,0,0.08,2,0,Dt,ParSim.vecSize,ParSim.LarmorFreq);
% SimPhantom.FIDs.NAA = Simulate_FID_Spectra(2.2154,4.65,0,0,0.08,3,0,Dt,ParSim.vecSize,ParSim.LarmorFreq);
SimPhantom.FIDs.NAA = Simulate_FID_Spectra(2.01,4.65,0,0,0.08,3,0,Dt,ParSim.vecSize,ParSim.LarmorFreq);
SimPhantom.FIDs.Cr = Simulate_FID_Spectra(3.02,4.65,0,0,0.08,3,0,Dt,ParSim.vecSize,ParSim.LarmorFreq);
SimPhantom.FIDs.Cho = Simulate_FID_Spectra(3.19,4.65,0,0,0.08,9,0,Dt,ParSim.vecSize,ParSim.LarmorFreq);

% % For Glutamate simulate a very simplified model with only one multiplet at 2.336 ppm. Tried to make it similar like in our 
% % Basis set for the FID sequence.
% dummy1 = Simulate_FID_Spectra(2.336,4.65,0,0,0.08,1,0,Dt,ParSim.vecSize,ParSim.LarmorFreq);
% dummy2 = Simulate_FID_Spectra(2.406,4.65,0,0,0.08,3/12,0,Dt,ParSim.vecSize,ParSim.LarmorFreq);
% dummy3 = Simulate_FID_Spectra(2.265,4.65,0,0,0.08,6/12,0,Dt,ParSim.vecSize,ParSim.LarmorFreq);
% SimPhantom.FIDs.Glu(2,:) = dummy1(2,:) + dummy2(2,:) + dummy3(2,:);
% clear dummy* Dt;

% Take only the FID
for CurField = transpose(fieldnames(SimPhantom.FIDs))
    SimPhantom.FIDs.(CurField{1}) = SimPhantom.FIDs.(CurField{1})(2,:);
end




%% Create Original Cartesian Data

if(prod(ParSim.MatSize)*2*8 / 2^30 < 10)                  % Data is < 10GB
    SimPhantom.csi = zeros(ParSim.MatSize);
    for CurField = transpose(fieldnames(SimPhantom.FIDs))
        SimPhantom.csi = SimPhantom.csi + myrepmat(SimPhantom.FIDs.(CurField{1}),ParSim.MatSize) .* ...
        myrepmat(SimPhantom.MetaboMaps.(CurField{1}),ParSim.MatSize);
    end
else
    
    error('\nCannot simulate data at once. Need to partition data! Please implement that in code!')
end

% SimPhantom.csi = ones(ParSim.MatSize);





