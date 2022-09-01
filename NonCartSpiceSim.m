% Description & Log:
% Changes vs previous version (SpiralSpiceSim_v4_4_TrajDeviations_B0Again3.m):
% * When resynthesizing the B0-expensive-Reco, use Reco_B0Corr.Phi instead of Reco.Phi! This made the reco working!
% * Bugfix: Fixed the bug that no noise was added to the D2-data.
% tbd: * Try if reco improves if I calculate Reco_B0Corr.Phi from a GroundTruth of size 64*64x400*159 instead of 64*64x304*159
% but then define the samplingOperator accordingly. IS THIS REALLY NECESSARY?
%
% Changes of SpiralSpiceSim_v4_4_TrajDeviations_B0Again3 vs SpiralSpiceSim_v4_4_TrajDeviations_B0Again2.m:
% In B0Again2 I additionally do an "expensive" Spice-reconstruction using rho of size 64*64x318*304. 
% Problem with B0Again2: It seems to work kind of, but the scaling of the "expensive" reco is wrong, and some vials in some 
% metabolites are very different...
% In B0Again2 I estimated the additional 304 time-points either by interpolation, or by using the ground truth
% Tried to simplify B0Again2 in 2 respects:
% * I have temporal interleaves in B0Again2, and I am not 100 % sure if I take care of it correctly. Therefore, changed to
% simulate no temporal interleaving for now.
% * Double-check if the Phi is correct that I use for the "expensive" SPICE-Reco




%% Housekeeping

clearvars
close all;

addpath('./dependencies');

quiet_flag = false;
TiltedTraj_flag = true;
B0_flag = true;
Correct4B0_flag = false;
ExpensiveB0Reco_flag = false;
FakeTilted_flag = false;
EqualMeasTime_flag = true;          % Change noise in conventional Spiral reco to simulate having same measurement time as SPICE


%% Constants & Parameters


SpSpice.SimPar.fov_overgrid = 1;
SpSpice.SimPar.fov_overgrid_SpiralSpice = 1;
SpSpice.SimPar.SpiralkSpaceSNR = 150;  % SNR = max(S(k,t=0))/(2*std(Noise). Use 0 for SNR of infinity (i.e. no noise added).
% SpSpice.SimPar.TempIntArtVec = [1 0.97*exp(1i*pi/16)];  % The Temporal Interleaves are multiplied with this, i.e. if this vector is
                                                      % [1 A*exp(1i*phi)], then the first TI is untouched, and the second time points
                                                      % are all multiplied by A, and phased with phase phi.
                                                      % BE CAREFUL THAT numel(SpSpice.SimPar.TempIntArtVec) = SpSpice.D2.Par.nTempInts
SpSpice.SimPar.TempIntArtVec = [1]; % [1 1] % (bstr): TI-CHANGE  % Simulate having no TI artifact
SpSpice.SimPar.Par.GradDelay_x = 3;
SpSpice.SimPar.Par.GradDelay_y = 4;



SpSpice.GroundTruth.Par.LarmorFreq = 123.223*10^6;
SpSpice.GroundTruth.Par.GammaBar = 4.2576 * 10^7;        % in Hz/T

% Cartesian Original Measurement Parameters
SpSpice.GroundTruth.Par.vecSize = 318; %318     % (bstr): TI-CHANGE
% SpSpice.GroundTruth.Par.vecSize = 1; % IMAGING-HACK: SIMULATE ONLY IMAGE, NO SPECTRO!

SpSpice.GroundTruth.Par.SpecDwellTimes = 1000000; % 1000000 % (bstr): TI-CHANGE
SpSpice.GroundTruth.Par.AcqTimePerTR = SpSpice.GroundTruth.Par.vecSize * SpSpice.GroundTruth.Par.SpecDwellTimes / 10^9; % seconds, 
% although the ADC is 0.32 s, the last two vecSize points are not measured
% It's the same also when having temporal interleaves, because then the Spec-dt doubles and the vecSize halves!



% SpSpice.GroundTruth.Par.SpecDwellTimes = SpSpice.GroundTruth.Par.AcqTimePerTR / SpSpice.GroundTruth.Par.vecSize * 10^9;            % in ns
SpSpice.GroundTruth.Par.ADC_Dt = SpSpice.GroundTruth.Par.SpecDwellTimes*1000;
SpSpice.GroundTruth.Par.MatSize = [64 64 1 SpSpice.GroundTruth.Par.vecSize];



% Define some Spiral Parameters
SpSpice.D2.Par.ADC_OverSamp = 2;
SpSpice.D2.Par.ADC_Dt = 10 / SpSpice.D2.Par.ADC_OverSamp;
SpSpice.D2.Par.nAngInts = 32;                               % Angular Interleaves
SpSpice.D2.Par.nTempInts = 1;    % 2 % (bstr): TI-CHANGE                           % Temporal Interleaves
% SpSpice.D2.Par.nTempInts = 1;                               % IMAGING-HACK: SIMULATE ONLY IMAGE, NO SPECTRO!
SpSpice.D2.Par.TrajFile = './InputMeasAndLogData/CONCEPTTrajFile_64x64_1000Hz_32nAI_1TI/';
fil = dir(SpSpice.D2.Par.TrajFile);
fil = {fil.name};
fil(cellfun(@isempty,regexp(fil,'\.m'))) = [];
SpSpice.D2.Par.TrajFile = strcat(SpSpice.D2.Par.TrajFile,fil);

run(SpSpice.D2.Par.TrajFile{1});
SpSpice.D2.Par.nSampPerTraj = NumberOfLoopPoints*SpSpice.D2.Par.ADC_OverSamp;   % Number of samples per interleaf
                                                                % (how many ADC-points make up one spiral). Simulate to have double!
SpSpice.D2.Par.vecSize = SpSpice.GroundTruth.Par.vecSize;
SpSpice.D2.Par.nRewPtsPerInt = 0;
SpSpice.D2.Par.MatSize = [SpSpice.D2.Par.nSampPerTraj SpSpice.D2.Par.nAngInts SpSpice.GroundTruth.Par.MatSize(end)]; % Dont consider nTIs
                                                                                                        % or temporal undersmapling here!
SpSpice.D2.Par.TimeUndersamplFactor = 2;

% Parameters of GroundTruth for calculating spiral. We might need to have a fine time-reaster for simulating a "tilted" trajectory (taking into account that the
% points along the spiral are measured at different time points
% Also, we might want to use a overgridding factor, thus increasing the MatSize/FoV
SpSpice.GroundTruth_ForSpir.Par = SpSpice.GroundTruth.Par;
% We will need a very fine raster in time to simulate that each spiral point is from a different time point
if(TiltedTraj_flag)
    % Since our trajectory is not instantaneuous, it's actually a little longer than the instantaneous phase encoding would be
    % It's longer by the length of one trajectory
%     SpSpice.GroundTruth_ForSpir.Par.AcqTimePerTR = SpSpice.GroundTruth_ForSpir.Par.AcqTimePerTR + ...
%     (SpSpice.D2.Par.nSampPerTraj+SpSpice.D2.Par.nRewPtsPerInt)/SpSpice.D2.Par.nTempInts*SpSpice.D2.Par.ADC_Dt/10^6;
%     SpSpice.GroundTruth_ForSpir.Par.vecSize = round(SpSpice.GroundTruth_ForSpir.Par.AcqTimePerTR / SpSpice.D2.Par.ADC_Dt*10^6);
    
    SpSpice.GroundTruth_ForSpir.Par.vecSize = round(SpSpice.GroundTruth_ForSpir.Par.AcqTimePerTR / SpSpice.D2.Par.ADC_Dt*10^6);
    SpSpice.GroundTruth_ForSpir.Par.vecSize = SpSpice.GroundTruth_ForSpir.Par.vecSize + ...
    (SpSpice.D2.Par.nTempInts-1)*(SpSpice.D2.Par.nSampPerTraj+SpSpice.D2.Par.nRewPtsPerInt)/SpSpice.D2.Par.nTempInts;
    SpSpice.GroundTruth_ForSpir.Par.AcqTimePerTR = SpSpice.GroundTruth_ForSpir.Par.vecSize * SpSpice.D2.Par.ADC_Dt/10^6;
end
SpSpice.GroundTruth_ForSpir.Par.MatSize(1:2) = SpSpice.SimPar.fov_overgrid*SpSpice.GroundTruth.Par.MatSize(1:2);
SpSpice.GroundTruth_ForSpir.Par.MatSize(4) = SpSpice.GroundTruth_ForSpir.Par.vecSize;
SpSpice.GroundTruth_ForSpir.Par.SpecDwellTimes = ...
SpSpice.GroundTruth_ForSpir.Par.AcqTimePerTR / SpSpice.GroundTruth_ForSpir.Par.vecSize * 10^9;   % in ns

% Parameters of GroundTruth for Calculating SpiralSpice. 
SpSpice.GroundTruth_ForSpSpice.Par = SpSpice.GroundTruth.Par;
SpSpice.GroundTruth_ForSpSpice.Par.MatSize(1:2) = SpSpice.SimPar.fov_overgrid_SpiralSpice*SpSpice.GroundTruth.Par.MatSize(1:2);

% Cartesian Elliptical. Against this data we compare our spiral reconstructions, bc it covers approximately the same (circular) k-space
SpSpice.CartEllip.Par = SpSpice.GroundTruth.Par;



% Parameters of Spi2Cart, i.e. spiral transformed to Cartesian image
SpSpice.Spi2Cart.Par = SpSpice.GroundTruth.Par;

% Parameters for GroundTruth Spiral, which is not influenced by B0, TiltedTrajectory, and noise
SpSpice.GroundTruth_Spiral.Par = SpSpice.D2.Par;

% Parameters for GroundTruth Spi2Cart, which is not influenced by B0, TiltedTrajectory, and noise
SpSpice.GroundTruth_Spi2Cart.Par = SpSpice.Spi2Cart.Par;



% Parameters of Spir2Cart_SpiralSpice
SpSpice.Reco.Par = SpSpice.GroundTruth.Par;

SpSpice.Reco_B0Corr.Par = SpSpice.GroundTruth.Par;








%% SIMULATION OF SPIRAL DATA:

%% Illustration

clear TIPts
% vecSize = 159; nTI = 1; TotTrajLength = 400; RewLength = 96; ADC_dt = 1; TrajLength = TotTrajLength - RewLength;
vecSize = SpSpice.D2.Par.vecSize; nTI = SpSpice.D2.Par.nTempInts; 
TotTrajLength = SpSpice.D2.Par.nRewPtsPerInt + SpSpice.D2.Par.nSampPerTraj; RewLength = SpSpice.D2.Par.nRewPtsPerInt;
ADC_dt = 1; TrajLength = SpSpice.D2.Par.nSampPerTraj;

PlotExtraPts = 2;
t = (0:TotTrajLength*vecSize/nTI+PlotExtraPts*TotTrajLength)*ADC_dt; FID = exp((-2*pi*1i*6.4/max(t)-6.4/max(t))*t); 
TIScheme = 1:nTI; TIScheme = repmat(TIScheme,[1 vecSize/nTI]);

for CurTI = 1:nTI
    TIPts{CurTI}.t = [];
    TIPts{CurTI}.FID = [];   
end
for CurvecS = 1:vecSize
    CurTI = TIScheme(CurvecS);
    
    TIPts{CurTI}.t = cat(2,TIPts{CurTI}.t,t( (CurvecS-1)*TotTrajLength/nTI+1 : (CurvecS-1)*TotTrajLength/nTI+TrajLength )); 
    TIPts{CurTI}.FID = cat(2,TIPts{CurTI}.FID,FID( (CurvecS-1)*TotTrajLength/nTI+1 : (CurvecS-1)*TotTrajLength/nTI+TrajLength ));     
end

figure; plot(t,real(FID),'k'), hold on,
scatter(TIPts{1}.t, TIPts{1}.FID,90,'b','x')
if(numel(TIPts) > 1)
    scatter(TIPts{2}.t, TIPts{2}.FID,30,'r','filled')
end
if(numel(TIPts) > 2)
    scatter(TIPts{3}.t, TIPts{3}.FID,30,'g','d','filled')
end

hold off
legend('FID','TI1','TI2','TI3')

NoOfPts = vecSize/nTI*TotTrajLength + (nTI-1)*TotTrajLength/nTI;


clear CurvecS CurTI

%% Simulate Phantom Data

ParSim = SpSpice.GroundTruth.Par;
run ./SimPhantomData3
SpSpice.GroundTruth.csi = SimPhantom.csi;
SpSpice.GroundTruth.Masks = SimPhantom.Masks;
clear ParSim SimPhantom;

if(TiltedTraj_flag)
    

    ParSim = SpSpice.GroundTruth_ForSpir.Par;
    run ./SimPhantomData3
    SpSpice.GroundTruth_ForSpir.csi = SimPhantom.csi;
    SpSpice.GroundTruth_ForSpir.Masks = SimPhantom.Masks;
    clear ParSim SimPhantom;
    SpSpice.GroundTruth_ForSpir.Par.MatSize = ...
    [SpSpice.GroundTruth_ForSpir.Par.MatSize(1:3) SpSpice.D2.Par.nTempInts SpSpice.D2.Par.nSampPerTraj SpSpice.D2.Par.vecSize/SpSpice.D2.Par.nTempInts];

        
    % Remove the time points that are not needed & Simulate Temporal Interleaving
    dummy = zeros([SpSpice.GroundTruth_ForSpir.Par.MatSize(1:4) prod(SpSpice.GroundTruth_ForSpir.Par.MatSize(5:6))]);
    KeepPts = cat(2,true([1 SpSpice.D2.Par.nSampPerTraj]),false([1 SpSpice.D2.Par.nRewPtsPerInt]));
    KeepPts = repmat(KeepPts,[1 SpSpice.D2.Par.vecSize/SpSpice.D2.Par.nTempInts]);
    HalfPtsPerTraj = (SpSpice.D2.Par.nSampPerTraj+SpSpice.D2.Par.nRewPtsPerInt)/SpSpice.D2.Par.nTempInts;
    for CurTI = 1:SpSpice.D2.Par.nTempInts
        KeepPts2 = cat(2,false([1 (CurTI-1)*HalfPtsPerTraj]),...
        KeepPts,false([1 (SpSpice.D2.Par.nTempInts-CurTI)*HalfPtsPerTraj]));
        dummy(:,:,:,CurTI,:) = SpSpice.GroundTruth_ForSpir.csi(:,:,:,KeepPts2);
        SpSpice.D2.TimeVec(CurTI,:) = (find(KeepPts2)-1) * SpSpice.GroundTruth_ForSpir.Par.SpecDwellTimes/10^9;
    end
    
%     SpSpice.GroundTruth_ForSpirExpB0 = SpSpice.GroundTruth_ForSpir;

    SpSpice.GroundTruth_ForSpir.csi = reshape(dummy,SpSpice.GroundTruth_ForSpir.Par.MatSize);
    SpSpice.GroundTruth_ForSpir.csi = permute(SpSpice.GroundTruth_ForSpir.csi,[1 2 3 5 4 6]);
    SpSpice.GroundTruth_ForSpir.Par.MatSize([4 5]) =  SpSpice.GroundTruth_ForSpir.Par.MatSize([5 4]);
    clear dummy;
    clear KeepPts;
    
    % Overgridding no longer supported in this case due to memory issues - would cause too high memory usage
    % LEGACY CODE: DONT USE GroundTruth_ForSpSpice.csi ANYLONGER
%     SpSpice.GroundTruth_ForSpSpice.csi = SpSpice.GroundTruth.csi;

    % Fake tilted trajectory, just copy the normal data (untilted) to check if reco is ok
    if(FakeTilted_flag)
%         SpSpice.GroundTruth_ForSpir.Par.MatSize([4 5]) =  SpSpice.GroundTruth_ForSpir.Par.MatSize([5 4]);    
%         SpSpice.GroundTruth_ForSpir.Par.MatSize = [SpSpice.GroundTruth_ForSpir.Par.MatSize 1 1];  % IMAGING-HACK: SIMULATE ONLY IMAGE, NO SPECTRO!
        SpSpice.GroundTruth_ForSpir.csi = ...
        reshape(SpSpice.GroundTruth.csi,[SpSpice.GroundTruth_ForSpir.Par.MatSize(1:3) 1 SpSpice.GroundTruth_ForSpir.Par.MatSize(5:6)]);
        SpSpice.GroundTruth_ForSpir.csi = repmat(SpSpice.GroundTruth_ForSpir.csi,[1 1 1 SpSpice.GroundTruth_ForSpir.Par.MatSize(4) 1 1]);
    end
    
    

else
    % Define 2x-FoV GroundTruths for calculating spiral and for usage in SpiralSpice
    
    % For Spiral
    SpSpice.GroundTruth_ForSpir.csi = ZerofillOrCutkSpace(SpSpice.GroundTruth.csi,[SpSpice.GroundTruth.Par.MatSize(1:2)*SpSpice.SimPar.fov_overgrid SpSpice.GroundTruth.Par.MatSize(3) SpSpice.GroundTruth.Par.MatSize(4)],0);
    SpSpice.GroundTruth_ForSpir.Par.MatSize = [SpSpice.GroundTruth_ForSpir.Par.MatSize(1:2)*SpSpice.SimPar.fov_overgrid_SpiralSpice SpSpice.GroundTruth_ForSpir.Par.MatSize(3) SpSpice.D2.Par.nTempInts 1 SpSpice.GroundTruth_ForSpir.Par.MatSize(end)/SpSpice.D2.Par.nTempInts];
    
    % Simulate Temporal Interleaving
    dummy = zeros([SpSpice.GroundTruth_ForSpir.Par.MatSize(1:4) prod(SpSpice.GroundTruth_ForSpir.Par.MatSize(5:6))]);
    for CurTI = 1:SpSpice.D2.Par.nTempInts
        dummy(:,:,:,CurTI,:) = SpSpice.GroundTruth_ForSpir.csi(:,:,:,CurTI:SpSpice.D2.Par.nTempInts:end);
    end
    SpSpice.GroundTruth_ForSpir.csi = reshape(dummy, SpSpice.GroundTruth_ForSpir.Par.MatSize);
    SpSpice.GroundTruth_ForSpir.csi = permute(SpSpice.GroundTruth_ForSpir.csi,[1 2 3 5 4 6]);
    SpSpice.GroundTruth_ForSpir.Par.MatSize([4 5]) =  SpSpice.GroundTruth_ForSpir.Par.MatSize([5 4]);
    clear dummy;
    
    
%     % LEGACY CODE: DONT USE GroundTruth_ForSpSpice.csi ANYLONGER
%     % For SpiralSpice
%     SpSpice.GroundTruth_ForSpSpice.csi = ZerofillOrCutkSpace(SpSpice.GroundTruth.csi,[SpSpice.GroundTruth.Par.MatSize(1:2)*SpSpice.SimPar.fov_overgrid_SpiralSpice SpSpice.GroundTruth.Par.MatSize(3) SpSpice.GroundTruth.Par.MatSize(4)],0);
%     SpSpice.GroundTruth_ForSpSpice.Par.MatSize = [SpSpice.GroundTruth_ForSpSpice.Par.MatSize(1:2)*SpSpice.SimPar.fov_overgrid_SpiralSpice SpSpice.GroundTruth_ForSpSpice.Par.MatSize(3) SpSpice.D2.Par.nTempInts 1 SpSpice.GroundTruth_ForSpSpice.Par.MatSize(end)/SpSpice.D2.Par.nTempInts];
% 
%     % Simulate Temporal Interleaving
%     dummy = zeros([SpSpice.GroundTruth_ForSpSpice.Par.MatSize(1:4) prod(SpSpice.GroundTruth_ForSpSpice.Par.MatSize(5:6))]);
%     for CurTI = 1:SpSpice.D2.Par.nTempInts
%         dummy(:,:,:,CurTI,:) = SpSpice.GroundTruth_ForSpSpice.csi(:,:,:,CurTI:SpSpice.D2.Par.nTempInts:end);
%     end
%     SpSpice.GroundTruth_ForSpSpice.csi = reshape(dummy, SpSpice.GroundTruth_ForSpSpice.Par.MatSize);
%     SpSpice.GroundTruth_ForSpSpice.csi = permute(SpSpice.GroundTruth_ForSpSpice.csi,[1 2 3 5 4 6]);
%     SpSpice.GroundTruth_ForSpSpice.Par.MatSize([4 5]) =  SpSpice.GroundTruth_ForSpSpice.Par.MatSize([5 4]);
%     clear dummy;
    
end

% % For determining the FudgeFactor in the density compensation function...
% SpSpice.GroundTruth.csi = ones(size(SpSpice.GroundTruth.csi));
% SpSpice.GroundTruth_ForSpir.csi = ones(size(SpSpice.GroundTruth_ForSpir.csi));


%% Analyze SV of GroundTruth Data

if(~quiet_flag)
    test = reshape(SpSpice.GroundTruth.csi,[prod(SpSpice.GroundTruth.Par.MatSize(1:end-1)) SpSpice.GroundTruth.Par.MatSize(end)]);
    [U, Sigma, V] = svd(test,'econ');
    figure; plot(diag(Sigma)), title('SVs of Cartesian Origianl Data')      % As expected we have rank 5 bc of water, NAA, Cr, Cho and Glx
end




%% Define CartEllip Data
% Remark: I could also use as the gold standard the BigFoV (zerofilled in iSpace) by using the BigFoV Gold standard and using a bigger
% Elliptical filter (bigger by the overgrid-factor). This gives slightly different results, but I think it shouldn't really matter.

SpSpice.CartEllip.csi = EllipticalFilter(SpSpice.GroundTruth.csi,[1 2],[1 1 1 ceil(SpSpice.GroundTruth.Par.MatSize(1)/2-1)],0);



%% Define the "Trajectory" of the GroundTruth_BigFoV

FoV = 0.220*SpSpice.SimPar.fov_overgrid; % in m
DeltaGM = 10^9/(FoV*SpSpice.GroundTruth.Par.GammaBar);           % in mT/m * us

% Calculate a Grid for PhaseEncoding Steps
% [bla_x, bla_y] = find(EllipticalFilter(ones(SpSpice.GroundTruth.Par.MatSize(1)*SpSpice.SimPar.fov_overgrid),[1 2],[1 1 1 ceil(SpSpice.GroundTruth.Par.MatSize(1)*SpSpice.SimPar.fov_overgrid/2-1)],1));
[bla_x, bla_y] = find(ones(SpSpice.GroundTruth_ForSpir.Par.MatSize(1)));

SpSpice.GroundTruth_ForSpir.Traj.GM(1,1,:) = bla_x - floor(SpSpice.GroundTruth.Par.MatSize(1)*SpSpice.SimPar.fov_overgrid/2) - 1; 
SpSpice.GroundTruth_ForSpir.Traj.GM(2,1,:) = bla_y - floor(SpSpice.GroundTruth.Par.MatSize(1)*SpSpice.SimPar.fov_overgrid/2) - 1; 

SpSpice.GroundTruth_ForSpir.Traj.GM = SpSpice.GroundTruth_ForSpir.Traj.GM * DeltaGM;

% Scale Trajectory
% maxR = max(max(abs(nsamp(1,:))), max(abs(nsamp(2,:))));
% We need to use the same maxR as we would use for fully-sampled data, bc otherwise we assign the wrong k-space positions to our
% k-space data (since just some k-space data are missing). If we calculate the max(max-stuff, we act as if the remaining k-space
% positions would be at the edge of the k-space (0.5,0.5), which they are not!
SpSpice.GroundTruth_ForSpir.Traj.maxR = 3416.4*1.0;          
SpSpice.GroundTruth_ForSpir.Traj.GM = SpSpice.GroundTruth_ForSpir.Traj.GM/SpSpice.GroundTruth_ForSpir.Traj.maxR/2;




%% Define the Spiral Trajectory

for jj = 1:numel(SpSpice.D2.Par.TrajFile)

    run(SpSpice.D2.Par.TrajFile{jj});
    
    CurAI = regexpi(SpSpice.D2.Par.TrajFile{jj},'SpatInterleaf[0-9]+\.m','match');
    CurAI = regexp(CurAI{1},'[0-9]+','match');
    CurAI = str2double(CurAI)+1;
    


    % The points in our trajectory are GradRasterTime apart from each other.
    % Simulate ADC_dt = GradientRaterTime/2, by replicating every gradient point. Because we sample finer than the GradRastertime,
    % the gradient value remains the same in every second ADC-point that we sample.
    % CurTraj = repmat(dGradientValues{1},[SpSpice.D2.Par.ADC_OverSamp 1]);
    % CurTraj = reshape(CurTraj,[1 numel(CurTraj)]);

    % Measured Trajectory
    GradDelay_x = SpSpice.SimPar.Par.GradDelay_x; GradDelay_y = SpSpice.SimPar.Par.GradDelay_y;
    CurTraj = dGradientValues{CurAI}; 
    CurTraj = CurTraj(1:NumberOfLoopPoints+NumberOfLaunTrackPoints);
    CurTraj = [zeros([1 2]) CurTraj CurTraj(NumberOfLaunTrackPoints+1:NumberOfLaunTrackPoints+3)];     % Append 2 0's --> max of 20 us gradient delay
    CurTraj = interp1(1:1:NumberOfLoopPoints+NumberOfLaunTrackPoints+5,CurTraj,1:0.1:NumberOfLoopPoints+NumberOfLaunTrackPoints+4.5); 
    CurTraj = complex(real(CurTraj((21-GradDelay_x):5:(end-GradDelay_x))), imag(CurTraj((21-GradDelay_y):5:(end-GradDelay_y))));
    % SpSpice.D2.TrajIdeal.GM = zeros([2 numel(CurTraj) SpSpice.D2.Par.nAngInts]);     % in mT/m * us
    blaa = -cumsum(CurTraj*dMaxGradAmpl*10/SpSpice.D2.Par.ADC_OverSamp);              % 5 is the ADC_dwelltime in us, GradientRaterTime = 2*ADC_dwelltime
    blaa = blaa((NumberOfLaunTrackPoints+2)*SpSpice.D2.Par.ADC_OverSamp+1:end);
    SpSpice.D2.TrajMeas.GM(:,:,CurAI) = [real(blaa); imag(blaa)];
%     SpSpice.D2.TrajMeas.GV(:,:,CurAI) = [real(CurTraj); imag(CurTraj)];

    
    
    % Ideal Trajectory
    % CurTraj = [0 dGradientValues{1}];
    % CurTraj = interp1(0:152,CurTraj,0.5:0.5:152.0);
    % New interpolation method: The last point has a problem, and I think it's not ok to append a 0 at beginning. So extrapolate
    CurTraj = dGradientValues{CurAI}; 
    CurTraj = CurTraj(1:NumberOfLoopPoints+NumberOfLaunTrackPoints);
    CurTraj = [CurTraj CurTraj(NumberOfLaunTrackPoints+1:NumberOfLaunTrackPoints+3)];     % Append 2 0's --> max of 20 us gradient delay    
    CurTraj = interp1(1:1:NumberOfLoopPoints+NumberOfLaunTrackPoints+3,CurTraj,1:0.5:NumberOfLoopPoints+NumberOfLaunTrackPoints+2.5); 
    % SpSpice.D2.TrajIdeal.GM = zeros([2 numel(CurTraj) SpSpice.D2.Par.nAngInts]);     % in mT/m * us
    blaa = -cumsum(CurTraj*dMaxGradAmpl*10/SpSpice.D2.Par.ADC_OverSamp);              % 5 is the ADC_dwelltime in us, GradientRaterTime = 2*ADC_dwelltime
    blaa = blaa((NumberOfLaunTrackPoints+2)*SpSpice.D2.Par.ADC_OverSamp+1:end);
    SpSpice.D2.TrajIdeal.GM(:,:,CurAI) = [real(blaa); imag(blaa)];
%     SpSpice.D2.TrajIdeal.GV(:,:,CurAI) = [real(CurTraj); imag(CurTraj)];












end

% Rotate spiral trajectory
Trajs = {'TrajIdeal','TrajMeas'};
for i=1:numel(Trajs)

    % Is it better to use the k-space radius of the Cartesian Fully sampled data?
    % Actually, it doesn't really matter, it will just simulate two different trajectories. Since we simulate the spiral data
    % and don't have them a priori (e.g. measured), we can use whatever trajectorie we want. However, of course if they are very
    % different than the Cartesian "trajectory", we will get high errors (e.g. if the spiral k-space extent is small, we can easily
    % grid to the spiral, but going back to the Cartesian is not fully possible, because we lost the data from the k-space periphery.
    % For now, I grid to the spiral that was simulated by poet, except that I use only half of the real ADC-dt. 
    % Therefore, I need to use the max k-space that we measure with a. fully sampled Cartesian k-space
    % That means that the 68x68 spiral is about the same size as the elliptical encoded 64x64

    SpSpice.D2.(Trajs{i}).maxR = 3416.4;
    % SpSpice.D2.(Trajs{i}).maxR = max(sqrt(SpSpice.D2.(Trajs{i}).GM(1,:).^2 + SpSpice.D2.(Trajs{i}).GM(2,:).^2))*1.0;        
    SpSpice.D2.(Trajs{i}).GM = SpSpice.D2.(Trajs{i}).GM/SpSpice.D2.(Trajs{i}).maxR/2;


end
clear CurTraj





%% Define the "Trajectory" of the CartEllip

FoV = 0.220; % in m
DeltaGM = 10^9/(FoV*SpSpice.CartEllip.Par.GammaBar);           % in mT/m * us

% Calculate a Grid for PhaseEncoding Steps
[bla_x, bla_y] = find(EllipticalFilter(ones(SpSpice.CartEllip.Par.MatSize(1)),[1 2],[1 1 1 ceil(SpSpice.CartEllip.Par.MatSize(1)/2-1)],1));
% [bla_x, bla_y] = find(ones(SpSpice.CartEllip.Par.MatSize(1)));

SpSpice.CartEllip.TrajIdeal.GM(1,1,:) = bla_x - floor(SpSpice.CartEllip.Par.MatSize(1)/2) - 1; 
SpSpice.CartEllip.TrajIdeal.GM(2,1,:) = bla_y - floor(SpSpice.CartEllip.Par.MatSize(1)/2) - 1; 

SpSpice.CartEllip.TrajIdeal.GM = SpSpice.CartEllip.TrajIdeal.GM * DeltaGM;

% Scale Trajectory
% maxR = max(max(abs(nsamp(1,:))), max(abs(nsamp(2,:))));
% We need to use the same maxR as we would use for fully-sampled data, bc otherwise we assign the wrong k-space positions to our
% k-space data (since just some k-space data are missing). If we calculate the max(max-stuff, we act as if the remaining k-space
% positions would be at the edge of the k-space (0.5,0.5), which they are not!
SpSpice.CartEllip.TrajIdeal.maxR = 3416.4*1.0;          
SpSpice.CartEllip.TrajIdeal.GM = SpSpice.CartEllip.TrajIdeal.GM/SpSpice.CartEllip.TrajIdeal.maxR/2;



%% Define the "Trajectory" of the GroundTruth_BigFoV_SpiralSpice

FoV = 0.220*SpSpice.SimPar.fov_overgrid_SpiralSpice; % in m
DeltaGM = 10^9/(FoV*SpSpice.GroundTruth.Par.GammaBar);           % in mT/m * us

% Calculate a Grid for PhaseEncoding Steps
[bla_x, bla_y] = find(ones(SpSpice.GroundTruth_ForSpSpice.Par.MatSize(1)));

SpSpice.GroundTruth_ForSpSpice.TrajIdeal.GM(1,1,:) = bla_x - floor(SpSpice.GroundTruth_ForSpSpice.Par.MatSize(1)/2) - 1; 
SpSpice.GroundTruth_ForSpSpice.TrajIdeal.GM(2,1,:) = bla_y - floor(SpSpice.GroundTruth_ForSpSpice.Par.MatSize(1)/2) - 1; 

SpSpice.GroundTruth_ForSpSpice.TrajIdeal.GM = SpSpice.GroundTruth_ForSpSpice.TrajIdeal.GM * DeltaGM;

% Scale Trajectory
% maxR = max(max(abs(nsamp(1,:))), max(abs(nsamp(2,:))));
% We need to use the same maxR as we would use for fully-sampled data, bc otherwise we assign the wrong k-space positions to our
% k-space data (since just some k-space data are missing). If we calculate the max(max-stuff, we act as if the remaining k-space
% positions would be at the edge of the k-space (0.5,0.5), which they are not!
SpSpice.GroundTruth_ForSpSpice.TrajIdeal.maxR = 3416.4*1.0;          
SpSpice.GroundTruth_ForSpSpice.TrajIdeal.GM = SpSpice.GroundTruth_ForSpSpice.TrajIdeal.GM/SpSpice.GroundTruth_ForSpSpice.TrajIdeal.maxR/2;





%% Plot Spiral and GroundTruth Trajectories

if(~quiet_flag)
    figure;
    scatter(squeeze(SpSpice.CartEllip.TrajIdeal.GM(1,1,:)),squeeze(SpSpice.CartEllip.TrajIdeal.GM(2,1,:)),'b'), hold on   
    scatter(squeeze(SpSpice.GroundTruth_ForSpir.Traj.GM(1,1,:)),squeeze(SpSpice.GroundTruth_ForSpir.Traj.GM(2,1,:)),'.k')
    for AngIntNo = 1:SpSpice.D2.Par.nAngInts
        scatter(squeeze(SpSpice.D2.TrajIdeal.GM(1,:,AngIntNo)), squeeze(SpSpice.D2.TrajIdeal.GM(2,:,AngIntNo)),'r')
        plot(squeeze(SpSpice.D2.TrajIdeal.GM(1,:,AngIntNo)), squeeze(SpSpice.D2.TrajIdeal.GM(2,:,AngIntNo)),'r')
        scatter(squeeze(SpSpice.D2.TrajMeas.GM(1,:,AngIntNo)), squeeze(SpSpice.D2.TrajMeas.GM(2,:,AngIntNo)),'g')
        plot(squeeze(SpSpice.D2.TrajMeas.GM(1,:,AngIntNo)), squeeze(SpSpice.D2.TrajMeas.GM(2,:,AngIntNo)),'g')
    end
    hold off
end






%% Define GroundTruth


nsamp = reshape(SpSpice.D2.TrajIdeal.GM,[2 SpSpice.D2.Par.nSampPerTraj*SpSpice.D2.Par.nAngInts]);
nsamp_Meas = reshape(SpSpice.D2.TrajMeas.GM,[2 SpSpice.D2.Par.nSampPerTraj*SpSpice.D2.Par.nAngInts]);

% sft Operator for ideal trajectory
sft2_Operator_ForSpir_TrajIdeal = sft2_Operator(...
transpose(squeeze(SpSpice.GroundTruth_ForSpir.Traj.GM)*size(SpSpice.GroundTruth_ForSpir.csi,1)),transpose(nsamp),1);

% sft Operator for measured trajectory. Used for creating the spiral data
sft2_Operator_ForSpir_TrajMeas = sft2_Operator(...
transpose(squeeze(SpSpice.GroundTruth_ForSpir.Traj.GM)*size(SpSpice.GroundTruth_ForSpir.csi,1)),transpose(nsamp_Meas),1);

% Define GroundTruth: Only take the first time point along the trajectory, having no additional phase
% and calculate the spiral data as if in the Non-Tilted case
SpSpice.GroundTruth_Spiral.csi_k = sft2_Operator_ForSpir_TrajIdeal * reshape(SpSpice.GroundTruth_ForSpir.csi(:,:,:,1,:,:), ...
[prod(SpSpice.GroundTruth_ForSpir.Par.MatSize(1:3)) prod(SpSpice.GroundTruth_ForSpir.Par.MatSize(5:end))]);
SpSpice.GroundTruth_Spiral.csi_k = reshape(SpSpice.GroundTruth_Spiral.csi_k,SpSpice.GroundTruth_Spiral.Par.MatSize);




%% Apply B0map (shift map)


SpSpice.GroundTruth_ForSpir.csi_NoB0 = SpSpice.GroundTruth_ForSpir.csi;
if(B0_flag)
    
    % Load B0 map
    load('./InputMeasAndLogData/B0Map.mat')            % In units of Hz
%     B0Map(abs(B0Map) > 50) = 0;
    
    % Make B0 map smooth again
    mask = imresize(SpSpice.GroundTruth.Masks.Vial{1},size(B0Map), 'nearest');
    VialMask = imresize(SpSpice.GroundTruth.Masks.VialMask,size(B0Map),'nearest');
    VialMask = MaskShrinkOrGrow(VialMask,6);
    mask_small = MaskShrinkOrGrow(mask,10,0);
    B0Map(mask==0 | (mask_small & VialMask & abs(B0Map) > 10)) = NaN;
    B0Map(B0Map > 60) = NaN;
    B0Map = inpaint_nans(B0Map).* mask;
    
    h = fspecial('gaussian',5,1.5);
    B0Map = imfilter(B0Map,h) .* mask;
    clear h;
%     figure; imagesc(B0Map)
    
    
    

    
    
    
    
%     % Another B0-Approach: Make an artifical B0-Map
%     B0Map_x = linspace(-123,123,64);
%     B0Map_y = linspace(123,-123,64);
%     [B0Map_x, B0Map_y] = ndgrid(B0Map_x, B0Map_y);
%     B0Map = (B0Map_x + B0Map_y);
%     B0Map = B0Map .* SpSpice.GroundTruth.Masks.Vial{1};
%     clear B0Map_*
    
    
    % Make B0-inhomogeneity stronger
    B0Map = 5*B0Map;
        
    % Apply to GroundTruth_ForSpir
    B0Map = imresize(B0Map,[size(SpSpice.GroundTruth_ForSpir.csi,1) size(SpSpice.GroundTruth_ForSpir.csi,2)]);
%     B0Map = B0Map .* SpSpice.GroundTruth.Masks.Vial{1};
    Sizey = size(SpSpice.GroundTruth_ForSpir.csi);
    
%     Sizey = [Sizey 1 1];% IMAGING-HACK: SIMULATE ONLY IMAGE, NO SPECTRO!
    Timey = reshape(SpSpice.D2.TimeVec,[size(SpSpice.GroundTruth_ForSpir.csi,5) Sizey(4) Sizey(6)]);
    SpSpice.D2.B0Mat = exp(myrepmat(-2*pi*1i*B0Map,Sizey) .* myrepmat(squeeze(Timey),Sizey));
    SpSpice.GroundTruth_ForSpir.csi = SpSpice.GroundTruth_ForSpir.csi .* SpSpice.D2.B0Mat;

    
    SpSpice.D2.B0Map = B0Map; clear B0Map B0Mat
else
    
    SpSpice.D2.B0Mat = ones(size(SpSpice.GroundTruth_ForSpir.csi));
    
end





% B0map_PointsToShift = randi([-B0map_Scale B0map_Scale],[MatSize(1) MatSize(2)]);
% ReferenceData = reshape(ReferenceData,[size(ReferenceData,1) size(ReferenceData,2) 1 size(ReferenceData,3)]);
% ReferenceData = squeeze(FrequencyAlignment(ReferenceData,B0map_PointsToShift,4,1));
% 
% % Simulate noisy B0map
% PointsToShift = -B0map_PointsToShift + randi([-B0map_Noise B0map_Noise],size(B0map_PointsToShift));
% 
% % Add some noise
% Noise = NoiseScale*randn(MatSize) + 1i*NoiseScale*randn(MatSize);
% PreSpice.NoB0Corr.InData = ReferenceData + Noise;









%% Define D1 Data

% For now, just use CartEllip data
SpSpice.D1 = SpSpice.CartEllip;


%% Create D2 Data: 
% Direct Fourier transform of the data to spiral using 2d slow Fourier transform (sft2)

% Probably we should use the SmallFoV, but it doesn't matter, because the BigFoV is just zerofilled (in iSpace), and these zeros
% do basically nothing!



% USE TrajMeas, NOT TrasIdeal !!!

% For DEBUGGING OF B0-CORR:
% SpSpice.D2.csi_k_Full = sft2_Operator_ForSpir_TrajMeas * reshape(SpSpice.GroundTruth_ForSpir.csi,[prod(SpSpice.GroundTruth_ForSpir.Par.MatSize(1:2)) prod(SpSpice.GroundTruth_ForSpir.Par.MatSize(3:end))]);
if(TiltedTraj_flag)
    % Apply sFT Operator, (spiral) point by point
    sft2_Operator_ForSpir_TrajMeas = reshape(sft2_Operator_ForSpir_TrajMeas,[SpSpice.D2.Par.nSampPerTraj SpSpice.D2.Par.nAngInts size(SpSpice.GroundTruth_ForSpir.Traj.GM,3)]);
    SpSpice.D2.csi_k = zeros([SpSpice.D2.Par.nSampPerTraj SpSpice.D2.Par.nAngInts SpSpice.D2.Par.vecSize]);
    for i = 1:size(sft2_Operator_ForSpir_TrajMeas,1)
        
        SpSpice.D2.csi_k(i,:,:) = squeeze(sft2_Operator_ForSpir_TrajMeas(i,:,:)) * ...
        reshape(SpSpice.GroundTruth_ForSpir.csi(:,:,:,i,:,:), [size(SpSpice.GroundTruth_ForSpir.csi,1)*... 
        size(SpSpice.GroundTruth_ForSpir.csi,2) prod(SpSpice.GroundTruth.Par.MatSize(3:4))]);
     end

    sft2_Operator_ForSpir_TrajMeas = reshape(sft2_Operator_ForSpir_TrajMeas,[SpSpice.D2.Par.nSampPerTraj*SpSpice.D2.Par.nAngInts size(SpSpice.GroundTruth_ForSpir.Traj.GM,3)]);
else   
    SpSpice.D2.csi_k = sft2_Operator_ForSpir_TrajMeas * reshape(SpSpice.GroundTruth_ForSpir.csi(:,:,:,1,:,:), ...
    [prod(SpSpice.GroundTruth_ForSpir.Par.MatSize(1:3)) prod(SpSpice.GroundTruth_ForSpir.Par.MatSize(5:end))]);
    SpSpice.D2.csi_k = reshape(SpSpice.D2.csi_k,SpSpice.GroundTruth_Spiral.Par.MatSize);
end

%%
% clearvars -except quiet_flag sft2_Operator_ForSpir_TrajIdeal sft2_Operator_ForSpir_TrajMeas SpSpice


%% Simulate Temporal Interleaving Artifact

for CurTI = 1:SpSpice.D2.Par.nTempInts
    SpSpice.D2.csi_k(:,:,CurTI:SpSpice.D2.Par.nTempInts:end) = ...
    SpSpice.D2.csi_k(:,:,CurTI:SpSpice.D2.Par.nTempInts:end)*SpSpice.SimPar.TempIntArtVec(CurTI);
end


%% Copy D2 to Spiral

SpSpice.Spiral = SpSpice.D2;


%% Add Noise to D2 Data


Signal = max(max(abs(SpSpice.D2.csi_k(:,:,1))));

NoiseStd = Signal / (2*SpSpice.SimPar.SpiralkSpaceSNR);        % from SpiralkSpaceSNR = Signal / (2*NoiseStd)
NoiseStd(isinf(NoiseStd)) = 0;

% % Comment in the following code for comparison, if it's important to always have the same IDENTICAL noise!
% if(exist('../InputMeasAndLogData/Noise.mat','file'))
%     load('../InputMeasAndLogData/Noise.mat')
% else
%     save('../InputMeasAndLogData/Noise.mat')
% end
% if(~all(size(Noise) == size(SpSpice.D2.csi_k))) 
%     Noise = NoiseStd*(randn(size(SpSpice.D2.csi_k)) + 1i*randn(size(SpSpice.D2.csi_k)));
% end
Noise = NoiseStd*(randn(size(SpSpice.D2.csi_k)) + 1i*randn(size(SpSpice.D2.csi_k)));

if(EqualMeasTime_flag)                                                      % If we measure x times longer, we have x times more signal
    Noise_Spiral = Noise * sqrt(0.9*SpSpice.D2.Par.TimeUndersamplFactor); % and sqrt(x) times more noise. Dividing by x to make the
else                                                                        % signal the same --> 1/sqrt(x) times more noise and same
    Noise_Spiral = Noise;                                                   % signal. To simulate same meas time --> need to multiply
end                                                                         % noise by sqrt(x), to achieve same SNR. 0.9: Spend ~10% on D1
    
SpSpice.Spiral.csi_k = SpSpice.Spiral.csi_k + Noise_Spiral;

SpSpice.D2.csi_k = SpSpice.D2.csi_k + Noise;

if(~quiet_flag)
    figure; plot(real(reshape(SpSpice.D2.csi_k(:,:,1),[1 prod(SpSpice.D2.Par.MatSize(1:2))]))); hold on; 
    plot(real(reshape(SpSpice.D2.csi_k(:,:,1),[1 prod(SpSpice.D2.Par.MatSize(1:2))])),'r'); hold off; legend('Noiseless','Noisy')
end


%% Define SamplingOperator for D2

SpSpice.D2.SamplingOperator = zeros([SpSpice.D2.Par.MatSize(1:end-1) SpSpice.D1.Par.MatSize(end)]);
SpSpice.D2.SamplingOperator_B0Corr = zeros([SpSpice.D2.Par.MatSize(1:end-1) SpSpice.D1.Par.MatSize(end)]);


% Randomly undersample in kx and ky
% test = zeros([80 80]);
% test(1:40,1:80) = 1;
% test_access = randperm(80*80);
% test2 = repmat(reshape(test(test_access),[80 80]),[1 1 1 SpSpice.D1.Par.MatSize(end)]);
% SpSpice.D2.SamplingOperator(logical(test2)) = 1;



SpSpice.D2.SamplingOperator(:,:,1:SpSpice.D2.Par.TimeUndersamplFactor:end) = 1;
SpSpice.D2.SamplingOperator_B0Corr(:,:,1:SpSpice.D2.Par.TimeUndersamplFactor:end) = 1;



%% SIM: Simulate Undersampling of D2, Create GoldStandard


SpSpice.D2.csi_k = SpSpice.D2.SamplingOperator .* SpSpice.D2.csi_k;



%% Define Results Size

SpSpice.Reco.Par.MatSize = SpSpice.GroundTruth.Par.MatSize;


























%% RECONSTRUCTION OF SPIRAL DATA (Gold Standard):


%% Try to correct for "tilted trajectory" by doing phaseroll

% Save the data for reconstructing the pseudo-pcg case
SpSpice.Spiral_NoTiltCorr = SpSpice.Spiral;
if(TiltedTraj_flag && ~FakeTilted_flag)

    nTI = SpSpice.D2.Par.nTempInts;
    vs = SpSpice.Spiral.Par.vecSize;
    ns = SpSpice.Spiral.Par.nSampPerTraj;
    nc = SpSpice.Spiral.Par.nAngInts;

    timeoffset = 0:(ns-1);
    timeoffset = repmat(transpose(timeoffset),[1 vs]);
    Freq = ((0:vs-1)/vs-0.5)/(SpSpice.Spiral.Par.nRewPtsPerInt + ns)*nTI;
    Freq = repmat(Freq,[ns 1]);    

    % This comes from:
    % timeoffset = (0:(ns-1))*SpSpice.Spiral.Par.ADC_Dt/10^6;
    % sBW = nTI/((SpSpice.Spiral.Par.nRewPtsPerInt + ns)*SpSpice.Spiral.Par.ADC_Dt/10^6);
    % Freq = -sBW/2 : sBW/vs : (sBW/2 - sBW/vs);
    % SpSpice.Spiral.Par.ADC_Dt/10^6 cancels out when calculating timeoffset * Freq and so can be omitted
    % the rest is basically the same (-sBW/2:sBW/vs:(sBW/2-sBW/vs) is equivalent to ((0:vs-1)/vs-0.5), and the other constants are
    % the same anyway
    
    phasecorr = exp(-2*1i*pi*timeoffset .* Freq);    
    phasecorr = repmat(reshape(phasecorr,[ns 1 vs]),[1 nc 1]);

    bla = SpSpice.Spiral.csi_k(:,:,1);
    SpSpice.Spiral.csi_k = fft(fftshift(conj(phasecorr).*fftshift(ifft(SpSpice.Spiral.csi_k,[],3),3),3),[],3);
    SpSpice.Spiral.csi_k(:,:,1) = bla;

    
    SpSpice.Grid.TiltTrajMat = reshape(phasecorr,[ns*nc vs]);   
    
else
    
	SpSpice.Grid.TiltTrajMat = ones([SpSpice.Spiral.Par.nSampPerTraj*SpSpice.Spiral.Par.nAngInts SpSpice.Spiral.Par.vecSize]);

end




%% Calculate B0-Correction of Spiral Data in Spatial Domain
if(B0_flag && TiltedTraj_flag)
    t   = (0:SpSpice.D2.Par.MatSize(1)-1)*SpSpice.D2.Par.ADC_Dt/10^6;
    t = repmat(t,[1 1 SpSpice.D2.Par.nAngInts]); t = t(:);
    SpSpice.Reco.B0CorrMat_Spatial = exp(myrepmat(-2*pi*1i*SpSpice.D2.B0Map(:),size(sft2_Operator_ForSpir_TrajIdeal)) .* myrepmat(t,size(sft2_Operator_ForSpir_TrajIdeal)));

    % SpSpice.Spi2Cart.B0CorrMat_Spatial = SpSpice.Reco.B0CorrMat_Spatial;
    % SpSpice.Spi2Cart_NoTiltCorr.B0CorrMat_Spatial = SpSpice.Reco.B0CorrMat_Spatial;

else
    SpSpice.Reco.B0CorrMat_Spatial = ones(size(sft2_Operator_ForSpir_TrajIdeal));
end
if(Correct4B0_flag)
    sft2_Operator_ForSpir_TrajIdeal_B0Corr = sft2_Operator_ForSpir_TrajIdeal .* SpSpice.Reco.B0CorrMat_Spatial;
else
    sft2_Operator_ForSpir_TrajIdeal_B0Corr = sft2_Operator_ForSpir_TrajIdeal;
end


%% Calculate Density Compensation According to Hoge1997 - Abrupt Changes

% FudgeFactor = 1.2743;     % For old trajectory
% FudgeFactor = 1.82;         % For new trajectory
% 
% v1 = SpSpice.Spiral.TrajIdeal.GM;
% SpSpice.Grid.SpirToGround.DCFPreG = zeros([1 size(v1,2) size(v1,3)]);
% for SpirPts = 2:size(v1,2)
%     SpSpice.Grid.SpirToGround.DCFPreG(1,SpirPts,:) = sqrt( v1(1,SpirPts,:).^2 + v1(2,SpirPts,:).^2 ) .* ...
%     abs( sqrt( v1(1,SpirPts,:).^2 + v1(2,SpirPts,:).^2 ) - sqrt( v1(1,SpirPts-1,:).^2 + v1(2,SpirPts-1,:).^2 ) );
% end
% SpSpice.Grid.SpirToGround.DCFPreG(isnan(SpSpice.Grid.SpirToGround.DCFPreG)) = 0;
% SpSpice.Grid.SpirToGround.DCFPreG = SpSpice.Grid.SpirToGround.DCFPreG(:) / max(SpSpice.Grid.SpirToGround.DCFPreG(:))*...
% 2*SpSpice.SimPar.fov_overgrid^2/FudgeFactor;  %
% % I dont know what these factors are. The 2*SpSpice.SimPar.fov_overgrid^2 I guessed. The FudgeFactor I got by inputting a image of ones
% % and seeing how it was scaled...

% clear v1


bla = (0.5:31.5) / 31.5;
bla = myrepmat(bla,[1 200 32]);
SpSpice.Grid.SpirToGround.DCFPreG = bla(:);

SpSpice.Grid.SpirToGround.DCFPreG_SpiralSpice = SpSpice.Grid.SpirToGround.DCFPreG;


%% Direct Fourier transform of the spiral data to Cartesian using 2d slow Fourier transform (sft2)




for InName = {'Spiral','Spiral_NoTiltCorr','GroundTruth_Spiral'}

    InData = SpSpice.(InName{:});
    OutName = regexprep(InName{:},'Spiral','Spi2Cart');
    
    % Noisy Data
    dummy = ...
    reshape(InData.csi_k,[InData.Par.nSampPerTraj*InData.Par.nAngInts InData.Par.vecSize]);   


    % Apply DCF
    SpSpice.Spi2Cart.csi_k = dummy .* myrepmat(SpSpice.Grid.SpirToGround.DCFPreG(:),size(dummy));
    clear dummy;

    if(~strcmpi(InName{:},'GroundTruth_Spiral'))       
        test2_dcfpre = sft2_Operator_ForSpir_TrajIdeal_B0Corr' * SpSpice.Spi2Cart.csi_k * size(SpSpice.GroundTruth_ForSpir.Traj.GM,3);
    else
        test2_dcfpre = sft2_Operator_ForSpir_TrajIdeal' * SpSpice.Spi2Cart.csi_k * size(SpSpice.GroundTruth_ForSpir.Traj.GM,3);        
    end

    dummy2 = reshape(test2_dcfpre,size(SpSpice.GroundTruth.csi) );            % NEEDS TO BE REPLACED FOR OVERGRIDDING!
    bla = (size(dummy2,1)-SpSpice.GroundTruth.Par.MatSize(1))/2+1;
    SpSpice.(OutName).csi = ...
    dummy2(bla:bla+SpSpice.GroundTruth.Par.MatSize(1)-1,bla:bla+SpSpice.GroundTruth.Par.MatSize(1)-1,:,:);

    clear test2_dcfpre bla dummy2

end


%% B0-Correction of Spiral Data in Spectral Domain (= Shifting Spectra)
if(B0_flag)
    t   = (0:SpSpice.Reco.Par.MatSize(end)-1)*SpSpice.Reco.Par.SpecDwellTimes/10^9;
    SpSpice.Reco.B0CorrMat_Spec = exp(myrepmat(2*pi*1i*SpSpice.D2.B0Map,SpSpice.Reco.Par.MatSize) .* myrepmat(t,SpSpice.Reco.Par.MatSize));


    % SpSpice.Spi2Cart.B0CorrMat_Spec = SpSpice.Reco.B0CorrMat_Spec;
    % SpSpice.Spi2Cart_NoTiltCorr.B0CorrMat_Spec = SpSpice.Reco.B0CorrMat_Spec;

    SpSpice.Spi2Cart_NoTiltCorr.csi = SpSpice.Spi2Cart_NoTiltCorr.csi .* SpSpice.Reco.B0CorrMat_Spec;
    SpSpice.Spi2Cart.csi = SpSpice.Spi2Cart.csi .* SpSpice.Reco.B0CorrMat_Spec;
else
    SpSpice.Reco.B0CorrMat_Spec = ones(SpSpice.Reco.Par.MatSize);
end




%% SPICE RECONSTRUCTION OF SPIRAL DATA:


%% Perform B0-Correction on D1, and Save B0-CorrectionMap for RecoData
% Obsolete for simulation

%% Remove NuisanceSignal of D1
% Obsolete for simulation


%% Extract Phi (Spectral Basis) from D1


CurData = SpSpice.D1.csi;       % To make life easier


% Rearrange Data to be  [ S(r1,t1)  ... S(r1,tN) ]
%                       [ S(r2,t1)  ... S(r2,tN) ]
%                                   ...
%                       [ S(rM,t1)  ... S(rM,tN) ]
SpSpice.D1.SVs = reshape(CurData,[prod(size_MultiDims(CurData,[1 2])) size(CurData,4)]); 

% Do an SVD
[SpSpice.D1.U, SpSpice.D1.SVs, SpSpice.D1.V] = svd(SpSpice.D1.SVs,'econ');


% % Estimate RankL Using Data + Noise
% % Obsolete, no noise in D1!
% % Gather noise from edge of image and end of FID
% SpiceSim.D1.NoiseHelpVars.NoiseVec = reshape(CurData,[1 size(CurData,1) size(CurData,2) 1 1 size(CurData,4)]);
% SpiceSim.D1.NoiseHelpVars.NoiseVec = GatherNoiseFromCSI(SpiceSim.D1.NoiseHelpVars.NoiseVec,1,200,200);
% 
% % Calculate Std
% SpiceSim.D1.NoiseHelpVars.NoiseVecStd = complex(std(real(SpiceSim.D1.NoiseHelpVars.NoiseVec)), std(imag(SpiceSim.D1.NoiseHelpVars.NoiseVec)));
% 
% % Create noise with same size as CurData and std as defined by SpiceSim.D1.NoiseHelpVars.NoiseVecStd
% SpiceSim.D1.NoiseHelpVars.Noise = complex(  real(SpiceSim.D1.NoiseHelpVars.NoiseVecStd)*randn(size(CurData)),imag(SpiceSim.D1.NoiseHelpVars.NoiseVecStd)*randn(size(CurData))  );
% 
% 
% % Do SVD of noise
% CurNoise = SpiceSim.D1.NoiseHelpVars.Noise .* SpiceSim.D1.tmask;
% 
% SpiceSim.D1.NoiseHelpVars.NoiseSVs = reshape( CurNoise,[prod(size_MultiDims(CurData,[1 2])) size(CurData,4)]);
% SpiceSim.D1.NoiseHelpVars.NoiseSVs = svd(SpiceSim.D1.NoiseHelpVars.NoiseSVs,'econ');
% 
% 
% % Plot Singular Values
% if(~quiet_flag)
%     figure; imagesc( ( std(real(CurData(:,:,:,end-300:end)),0,4) + std(imag(CurData(:,:,:,end-300:end)),0,4) ) / 2)
%     figure; plot(diag(SpiceSim.D1.SVs),'b'); hold on; plot(SpiceSim.D1.NoiseHelpVars.NoiseSVs,'r'); hold off;
%     title('Singular Values B0Corr'); legend('SVs(B0Corr)','SVs(Noise)')
% end
%
%
% clearvars -except SpiceSim quiet_flag
% RankL_Suggested = find(diag(SpiceSim.D1.SVs) < SpiceSim.D1.NoiseHelpVars.NoiseSVs(1),1) %#ok
% RankL = RankL_Suggested;
% 
% RankL_max = 30;
% if(RankL > RankL_max)
%     fprintf('\nRestricting rank from %d to %d',RankL, RankL_max);
%     RankL(RankL > RankL_max) = RankL_max;
% end

RankL = 4;
% RankL = 1;% IMAGING-HACK: SIMULATE ONLY IMAGE, NO SPECTRO!
SpSpice.Reco.Phi = (SpSpice.D1.V(:,1:RankL)');
SpSpice.Reco.RankL = RankL;




%% SIM: Create Denoised Spi2Cart Data for Comparison with SPICE-Reco

SpSpice.Spi2Cart_Denoised.csi = SpSpice.Spi2Cart.csi;
SpSpice.Spi2Cart_Denoised.Par = SpSpice.Spi2Cart.Par;

% Perfrom B0-correction and Nuisance Removal on GoldStandard
% Obsolete

% Nuisance Removal 
% Obsolete



% Truncated SVD Gold Standard
CurData = SpSpice.Spi2Cart_Denoised.csi;

% Apply mask
% Obsolete
SVs = CurData; clear CurData;
SVs = reshape(SVs,[prod(size_MultiDims(SVs,1:3)) size(SVs,4)]);

% Do an SVD
[U, SVs, V] = svd(SVs,'econ');
% SpiceSim.GoldStandard.U = U; SpiceSim.GoldStandard.SVs = SVs; SpiceSim.GoldStandard.V = V; 

RankL_GoldStandard = RankL;         % Chao said they should be the same for comparison reasons...

% Resynthesize Data with Truncated SVD
SpSpice.Spi2Cart_Denoised.csi = U(:,1:RankL_GoldStandard) * SVs(1:RankL_GoldStandard,1:RankL_GoldStandard) * (V(:,1:RankL_GoldStandard)');
SpSpice.Spi2Cart_Denoised.csi = reshape(SpSpice.Spi2Cart_Denoised.csi, SpSpice.Spi2Cart_Denoised.Par.MatSize);

clearvars -except SpSpice quiet_flag RankL nsamp* *_flag



%% Remove NuisanceSignal of D2
% Obsolete for simulation




%% Reconstruct Spatial Map UTS (= U Times Sigma)


tic_PCG = tic;

% Define Input Data
% The ratio of the FoV-overgrids is necessary, because the data is scaled by the matrix size when going Cartesian --> Spiral, but
% the Cartesian data is zerofilled in image domain, and those zeros dont contribute any signal...
CartSize = SpSpice.GroundTruth_ForSpSpice.Par.MatSize(1);
D2 = (SpSpice.SimPar.fov_overgrid/SpSpice.SimPar.fov_overgrid_SpiralSpice)^2*SpSpice.D2.csi_k;
sft2_Operator_ForSpSpice = sft2_Operator(transpose(squeeze(SpSpice.GroundTruth_ForSpSpice.TrajIdeal.GM)*CartSize),transpose(nsamp),1);

if(Correct4B0_flag)
    sft2_Operator_ForSpSpice = sft2_Operator_ForSpSpice .* SpSpice.Reco.B0CorrMat_Spatial;
end

% B0CorrMat = ones([prod(SpSpice.Reco.Par.MatSize(1:end-1)) SpSpice.Reco.Par.MatSize(end)]);
B0CorrMat = conj(reshape(SpSpice.Reco.B0CorrMat_Spec,[prod(SpSpice.Reco.Par.MatSize(1:end-1)) SpSpice.Reco.Par.MatSize(end)]));

SpSpice.D2.SamplingOperator = reshape(SpSpice.D2.SamplingOperator,[prod(SpSpice.D2.Par.MatSize(1:2)) SpSpice.D2.Par.MatSize(end)]);

D2 = reshape(D2, [numel(D2) 1]);
D2 = SpSpice_SynthesizeMeasData_v4_4_TrajDeviations(...
'Transj',D2,SpSpice.Reco.Phi,SpSpice.D2.SamplingOperator,B0CorrMat,SpSpice.Grid.SpirToGround.DCFPreG_SpiralSpice,sft2_Operator_ForSpSpice,SpSpice.Grid.TiltTrajMat );


funh = @(x) SpSpice_SynthesizeMeasData_v4_4_TrajDeviations('Transj', SpSpice_SynthesizeMeasData_v4_4_TrajDeviations('NoTransj',x,...
SpSpice.Reco.Phi,SpSpice.D2.SamplingOperator,B0CorrMat,SpSpice.Grid.SpirToGround.DCFPreG_SpiralSpice,sft2_Operator_ForSpSpice,SpSpice.Grid.TiltTrajMat),      ...
SpSpice.Reco.Phi,SpSpice.D2.SamplingOperator,B0CorrMat,SpSpice.Grid.SpirToGround.DCFPreG_SpiralSpice,sft2_Operator_ForSpSpice,SpSpice.Grid.TiltTrajMat);


% InitGuess = reshape(SpiceSim.ResynthD2.U,[numel(SpiceSim.ResynthD2.U) 1]);
% InitGuess = InitGuess + 0.04*max(abs(InitGuess(:)))*randn(size(InitGuess));
% InitGuess = ones([numel(SpiceSim.D1.U) 1]);
InitGuess = [];
% InitGuess = randn([numel(SpiceSim.D1.U) 1]);


fprintf('\n\nStart pcg . . .')
SpSpice.Reco.UTS = pcg(funh,D2,1E-6,20,[],[],InitGuess);
SpSpice.Reco.UTS = reshape(SpSpice.Reco.UTS,[SpSpice.Reco.Par.MatSize(1:end-1) RankL]);

SpSpice.D2.SamplingOperator = reshape(SpSpice.D2.SamplingOperator,SpSpice.D2.Par.MatSize);



% Resynthesize Data Using Estimated UTS
SpSpice.Reco.UTS = reshape(SpSpice.Reco.UTS,[prod(SpSpice.Reco.Par.MatSize(1:end-1)) RankL]);
SpSpice.Reco.csi = SpSpice.Reco.UTS(:,1:RankL) * SpSpice.Reco.Phi;
SpSpice.Reco.csi = reshape(SpSpice.Reco.csi,[SpSpice.Reco.Par.MatSize(1:end-1) SpSpice.D1.Par.MatSize(end)]);

% Cut k-space to measured k-space extent
% bla = FFTOfMRIData(SpSpice.Reco.csi,1,[1 2],1,1,0);
% Filter = EllipticalFilter(ones(size(bla)),[1 2],[1 1 1 ceil(SpSpice.GroundTruth.Par.MatSize(1)/2-1)],1);
% SpSpice.Reco.csi = FFTOfMRIData(bla.*Filter,1,[1 2],0,1,0);

fprintf('\ntook %f s',toc(tic_PCG))

clearvars -except SpSpice quiet_flag RankL nsamp *_flag



%% Reconstruct Spatial Map UTS (= U Times Sigma) Correcting for B0 Expensively!


tic_PCG_B0Corr = tic;

%% Define or Reshape Operators
sft2_Operator_ForSpSpice = sft2_Operator(transpose(squeeze(SpSpice.GroundTruth_ForSpSpice.TrajIdeal.GM)*SpSpice.GroundTruth_ForSpSpice.Par.MatSize(1)),transpose(nsamp),1);

SpSpice.D2.SamplingOperator_B0Corr = false([SpSpice.D2.Par.MatSize(1:end-1) SpSpice.D2.Par.MatSize(1)*SpSpice.D2.Par.MatSize(end)]);

for ii = 1:SpSpice.D2.Par.MatSize(1)
    SpSpice.D2.SamplingOperator_B0Corr(ii,:,ii:(SpSpice.D2.Par.MatSize(1)*SpSpice.D2.Par.TimeUndersamplFactor):ii+SpSpice.D2.Par.MatSize(1)*(SpSpice.D2.Par.MatSize(end)-1)) = true; 
end
SpSpice.D2.SamplingOperator_B0Corr = reshape(SpSpice.D2.SamplingOperator_B0Corr,[prod(SpSpice.D2.Par.MatSize(1:2)) size(SpSpice.D2.SamplingOperator_B0Corr,3)]);


%SpSpice.D2.B0Mat = SpSpice.D2.B0Mat(:,:,:,:,:,1:158);
SpSpice.D2.B0Mat = ...
reshape(SpSpice.D2.B0Mat,[prod(SpSpice.GroundTruth_ForSpir.Par.MatSize(1:2)) prod(size_MultiDims(SpSpice.D2.B0Mat,[4 5 6]))]);
% SpSpice.D2.B0Mat = SpSpice.D2.B0Mat(:,1:size(D2,3));



%% Get Phi with Fine Time Steps


if(ExpensiveB0Reco_flag)
    
    
    % Option1: Interpolate Phi
    % 
    % % sft2_Operator_ForSpSpice = sft2_Operator_ForSpSpice .* SpSpice.Reco.B0CorrMat_Spatial;
    % 
    % 
    % ZerfillFac = 200;
    % NoOfTimePts = 318; 
    % 
    % % By zerofilling in frequency domain
    % bla = SpSpice.Reco.Phi; 
    % 
    % % Mirrored Zerofilling
    % bla = cat(2,conj(bla(:,end:-1:2)),bla);
    % bla2 = zeros([size(bla,1) size(bla,2)+1]);
    % for ii=1:size(bla,1)
    %     bla2(ii,:) = interp1(1:numel(bla(ii,:)),bla(ii,:),0:numel(bla(ii,:)),'linear','extrap');
    % end
    % bla = fftshift(fft(fftshift(bla2*200,2),[],2),2); clear bla2
    % bla = ZerofillOrCutkSpace(bla,[4 2*ZerfillFac*NoOfTimePts],0);  % IT WAS: (NoOfTimePts-1). WHY -1?
    % bla = ifftshift(ifft(ifftshift(bla,2),[],2),2);
    % bla = bla(:,end/2+1:end); 
    % 
    % % % Normal Zerofilling
    % % bla = fftshift(fft(bla,[],2),2); clear bla2
    % % bla = ZerofillOrCutkSpace(bla,[4 ZerfillFac*(NoOfTimePts-1)],0); 
    % % bla = ifft(ifftshift(bla,2),[],2);
    % 
    % 
    % % bla = bla(:,1:315*200+304);
    % 
    % % Shuffle time points of Phi, because we measure some time points twice 
    % bla2 = zeros([4 316*304]);
    % for ii = 1:316
    %    bla2(:,(ii-1)*304+1:ii*304) = bla(:,(ii-1)*200+1:(ii-1)*200+304);
    % end
    % bla = bla2; clear bla2;
    % 
    % % for ii = 1:size(bla,1)
    % %     bla2(ii,:) = interp1(1:size(bla,2),bla(ii,:),1:size(bla,2)+200,'cubic','extrap');
    % % end
    % % bla = bla2;
    % 
    % SpSpice.Reco_B0Corr.Phi = bla;
    % 
    % % % Interpolate in time
    % % % Using spline
    % % bla = SpSpice.Reco.Phi;
    % % SpSpice.Reco_B0Corr.Phi = zeros([size(bla,1) ZerfillFac*NoOfTimePts]);
    % % for ii = 1:size(bla,1)
    % %     SpSpice.Reco_B0Corr.Phi(ii,1:ZerfillFac*(NoOfTimePts-1)) = interp1(1:(NoOfTimePts-1),bla(ii,:),1:1/ZerfillFac:(NoOfTimePts-1/ZerfillFac),'spline','extrap');
    % % end
    % 
    % 
    % clear bla bla2




    % Option2: Cheat and use GroundTruth to estimate correct Phi
    [UU SS VV] = svd(reshape(SpSpice.GroundTruth_ForSpir.csi_NoB0(:,:,:,:,:,1:end),[prod(SpSpice.GroundTruth_ForSpir.Par.MatSize(1:2)) prod(SpSpice.GroundTruth_ForSpir.Par.MatSize(3:end))]),'econ');
    figure; plot(diag(SS))
    UTS2 = UU(:,1:4) * SS(1:4,1:4);
    Phi2 = (VV(:,1:4))'; save('./UTS2_ConcRings.mat','UTS2','Phi2');
    load('./UTS2_ConcRings.mat')
    SpSpice.Reco_B0Corr.Phi = Phi2;


    % Option3: Perform an iterative reconstruction using D1-data to Get Phi

end


%% Memory Clearing

SpSpice.GroundTruth_ForSpir.csi = [];
SpSpice.GroundTruth_ForSpir.csi_NoB0 = [];



%% Get D2 with Fine Timesteps


if(ExpensiveB0Reco_flag)
    % Option1: Every point we haven't measured, set to 0 (zerofill D2)
    CartSize = SpSpice.GroundTruth_ForSpSpice.Par.MatSize(1);

    D2 = zeros([SpSpice.D2.Par.MatSize(1:2) size(SpSpice.Reco_B0Corr.Phi,2)]);
    for ii = 1:size(D2,1)
        D2(ii,:,ii:SpSpice.D2.Par.MatSize(1):ii+(SpSpice.D2.Par.MatSize(3)-1)*SpSpice.D2.Par.MatSize(1)) = SpSpice.D2.csi_k(ii,:,1:SpSpice.D2.Par.MatSize(3));
    end
    D2 = (SpSpice.SimPar.fov_overgrid/SpSpice.SimPar.fov_overgrid_SpiralSpice)^2*D2;




    % Option2: Interpolate D2

    % bla = D2;
    % % Mirrored Zerofilling
    % bla = cat(3,conj(bla(:,:,end:-1:2)),bla);
    % bla2 = zeros([size(bla,1) size(bla,2) size(bla,3)+1]);
    % for ii=1:size(bla,1)
    %     for jj=1:size(bla,2)
    %         bla2(ii,jj,:) = interp1(1:size(bla,3),squeeze(bla(ii,jj,:)),0:size(bla,3),'linear','extrap');
    %     end
    % end
    % bla = fftshift(fft(fftshift(bla2*200,3),[],3),3); clear bla2
    % bla = ZerofillOrCutkSpace(bla,[size(bla,1) size(bla,2) 2*ZerfillFac*(NoOfTimePts-1)],0); 
    % bla = ifftshift(ifft(ifftshift(bla,3),[],3),3);
    % bla = bla(:,:,end/2+1:end); 
    % 
    % % bla2 = zeros([size(bla,1) size(bla,2) size(bla,3)+200]);
    % % for ii = 1:size(bla,1)
    % % 	for jj=1:size(bla,2)
    % %         bla2(ii,jj,:) = interp1(1:size(bla,3),squeeze(bla(ii,jj,:)),1:size(bla,3)+200,'cubic','extrap');
    % % 	end
    % % end
    % % bla = bla2;
    % % size(bla)
    % 
    % D2 = reshape(bla, [numel(bla) 1]);
    % clear bla bla2

end



%% Perform Reconstruction


if(ExpensiveB0Reco_flag)

    % SampOp = reshape(repmat(reshape(SpSpice.D2.SamplingOperator_B0Corr,[4560 1 318]),[1 ZerfillFac 1]),[4560 ZerfillFac*318]);
    % SpSpice.D2.SamplingOperator_B0Corr = reshape(SpSpice.D2.SamplingOperator_B0Corr,[prod(SpSpice.D2.Par.MatSize(1:2)) SpSpice.D2.Par.MatSize(end)]);

    % reshape
    D2 = reshape(D2,[size(D2,1)*size(D2,2) size(D2,3)]);
    SpSpice.D2.SamplingOperator_B0Corr = reshape(SpSpice.D2.SamplingOperator_B0Corr,[prod(SpSpice.D2.Par.MatSize(1:2)) size(D2,2)]);



    % Estimate UTS
    D2 = SpSpice_SynthesizeMeasData_v4_4_TrajDeviations_Expensive2(...
    'Transj',D2,SpSpice.Reco_B0Corr.Phi,SpSpice.D2.SamplingOperator_B0Corr,SpSpice.D2.B0Mat,SpSpice.Grid.SpirToGround.DCFPreG_SpiralSpice,sft2_Operator_ForSpSpice,SpSpice.Grid.TiltTrajMat );

    funh = @(x) SpSpice_SynthesizeMeasData_v4_4_TrajDeviations_Expensive2('Transj', SpSpice_SynthesizeMeasData_v4_4_TrajDeviations_Expensive2('NoTransj',x,...
    SpSpice.Reco_B0Corr.Phi,SpSpice.D2.SamplingOperator_B0Corr,SpSpice.D2.B0Mat,SpSpice.Grid.SpirToGround.DCFPreG_SpiralSpice,sft2_Operator_ForSpSpice,SpSpice.Grid.TiltTrajMat),      ...
    SpSpice.Reco_B0Corr.Phi,SpSpice.D2.SamplingOperator_B0Corr,SpSpice.D2.B0Mat,SpSpice.Grid.SpirToGround.DCFPreG_SpiralSpice,sft2_Operator_ForSpSpice,SpSpice.Grid.TiltTrajMat);


    % InitGuess = reshape(SpiceSim.ResynthD2.U,[numel(SpiceSim.ResynthD2.U) 1]);
    % InitGuess = InitGuess + 0.04*max(abs(InitGuess(:)))*randn(size(InitGuess));
    % InitGuess = ones([numel(SpiceSim.D1.U) 1]);
    InitGuess = [];
    % InitGuess = randn([numel(SpiceSim.D1.U) 1]);


    fprintf('\n\nStart pcg . . .')
    SpSpice.Reco_B0Corr.UTS = pcg(funh,D2,1E-6,10,[],[],InitGuess);
    SpSpice.Reco_B0Corr.UTS = reshape(SpSpice.Reco_B0Corr.UTS,[SpSpice.Reco_B0Corr.Par.MatSize(1:end-1) RankL]);

    SpSpice.D2.SamplingOperator = reshape(SpSpice.D2.SamplingOperator,SpSpice.D2.Par.MatSize);



    % Resynthesize Data Using Estimated UTS
    SpSpice.Reco_B0Corr.UTS = reshape(SpSpice.Reco_B0Corr.UTS,[prod(SpSpice.Reco_B0Corr.Par.MatSize(1:end-1)) RankL]);
    SpSpice.Reco_B0Corr.csi = SpSpice.Reco_B0Corr.UTS(:,1:RankL) * SpSpice.Reco_B0Corr.Phi;
    SpSpice.Reco_B0Corr.csi = SpSpice.Reco_B0Corr.csi(:,1:SpSpice.D2.Par.MatSize(1):end);
    SpSpice.Reco_B0Corr.csi = reshape(SpSpice.Reco_B0Corr.csi,[SpSpice.Reco_B0Corr.Par.MatSize(1:end-1) SpSpice.D1.Par.MatSize(end)]);

    % Cut k-space to measured k-space extent
    % bla = FFTOfMRIData(SpSpice.Reco_B0Corr.csi,1,[1 2],1,1,0);
    % Filter = EllipticalFilter(ones(size(bla)),[1 2],[1 1 1 ceil(SpSpice.GroundTruth.Par.MatSize(1)/2-1)],1);
    % SpSpice.Reco_B0Corr.csi = FFTOfMRIData(bla.*Filter,1,[1 2],0,1,0);

    fprintf('\ntook %f s',toc(tic_PCG_B0Corr))

    clearvars -except SpSpice quiet_flag RankL nsamp


    % SpSpice.Reco = SpSpice.Reco_B0Corr;


else
    SpSpice.Reco_B0Corr = SpSpice.Reco;
    
end























%% DISPLAY RESULTS

% vec = 1;% IMAGING-HACK: SIMULATE ONLY IMAGE, NO SPECTRO!
vec = 2;
maxscale = 1.5*max(max(abs(SpSpice.GroundTruth_Spi2Cart.csi(:,:,:,vec))));
mask = SpSpice.GroundTruth.Masks.Vial{1};
tmask = repmat(mask,[1 1 1 SpSpice.GroundTruth_Spi2Cart.Par.MatSize(end)]);
% mask = SpSpice.GroundTruth.Masks.VialMask;
% mask = ones(size(SpSpice.GroundTruth.Masks.VialMask));


% GroundTruth_Spi2Cart vs GroundTruth
RelErr = ...
(abs(SpSpice.GroundTruth_Spi2Cart.csi - (SpSpice.GroundTruth.csi))) ;
RelErr(tmask == 0) = NaN;
RelError.SpirGroundTruth_vs_GroundTruth = nanmean_own(RelErr(:));
RelErrorMap.SpirGroundTruth_vs_GroundTruth = RelErr(:,:,:,vec);


% Spiral vs GroundTruth
RelErr = ...
(abs(SpSpice.Spi2Cart.csi - (SpSpice.GroundTruth.csi))) ;
RelErr(tmask == 0) = NaN;
RelError.Spiral_vs_GroundTruth = nanmean_own(RelErr(:));
RelErrorMap.Spiral_vs_GroundTruth = RelErr(:,:,:,vec);


% SpiralDenoised vs GroundTruth
RelErr = ...
(abs(SpSpice.Spi2Cart_Denoised.csi - (SpSpice.GroundTruth.csi))) ;
RelErr(tmask == 0) = NaN;
RelError.SpiralDenoised_vs_GroundTruth = nanmean_own(RelErr(:));
RelErrorMap.SpiralDenoised_vs_GroundTruth = RelErr(:,:,:,vec);


% SpiralSpice vs GroundTruth
RelErr = ...
(abs(SpSpice.Reco.csi - (SpSpice.GroundTruth.csi))) ;
RelErr(tmask == 0) = NaN;
RelError.SpiralSpice_vs_GroundTruth = nanmean_own(RelErr(:));
RelErrorMap.SpiralSpice_vs_GroundTruth = RelErr(:,:,:,vec);


% SpiralSpiceB0Corr vs GroundTruth
RelErr = ...
(abs(SpSpice.Reco_B0Corr.csi - (SpSpice.GroundTruth.csi))) ;
RelErr(tmask == 0) = NaN;
RelError.SpiralSpiceB0Corr_vs_GroundTruth = nanmean_own(RelErr(:));
RelErrorMap.SpiralSpiceB0Corr_vs_GroundTruth = RelErr(:,:,:,vec);





% Spiral vs GroundTruth_Spi2Cart
RelErr = ...
(abs(SpSpice.Spi2Cart.csi - (SpSpice.GroundTruth_Spi2Cart.csi))) ;
RelErr(tmask == 0) = NaN;
RelError.Spiral_vs_SpirGroundTruth = nanmean_own(RelErr(:));
RelErrorMap.Spiral_vs_SpirGroundTruth = RelErr(:,:,:,vec);


% SpiralDenoised vs GroundTruth_Spi2Cart
RelErr = ...
(abs(SpSpice.Spi2Cart_Denoised.csi - (SpSpice.GroundTruth_Spi2Cart.csi))) ;
RelErr(tmask == 0) = NaN;
RelError.SpiralDenoised_vs_SpirGroundTruth = nanmean_own(RelErr(:));
RelErrorMap.SpiralDenoised_vs_SpirGroundTruth = RelErr(:,:,:,vec);


% SpiralSpice vs GroundTruth_Spi2Cart
RelErr = ...
(abs(SpSpice.Reco.csi - (SpSpice.GroundTruth_Spi2Cart.csi))) ;
RelErr(tmask == 0) = NaN;
RelError.SpiralSpice_vs_SpirGroundTruth = nanmean_own(RelErr(:));
RelErrorMap.SpiralSpice_vs_SpirGroundTruth = RelErr(:,:,:,vec);


% SpiralSpiceB0Corr vs GroundTruth_Spi2Cart
RelErr = ...
(abs(SpSpice.Reco_B0Corr.csi - (SpSpice.GroundTruth_Spi2Cart.csi))) ;
RelErr(tmask == 0) = NaN;
RelError.SpiralSpiceB0Corr_vs_SpirGroundTruth = nanmean_own(RelErr(:));
RelErrorMap.SpiralSpiceB0Corr_vs_SpirGroundTruth = RelErr(:,:,:,vec);




% if(~quiet_flag)
    figure; 
    subplot(3,2,1); imagesc(RelErrorMap.SpirGroundTruth_vs_GroundTruth .* mask), colorbar, title('RelError SpiralGroundTruth vs GroundTruth')
    subplot(3,2,2); imagesc(RelErrorMap.Spiral_vs_SpirGroundTruth .* mask), colorbar, title(['RelError Spiral vs SpiralGroundTruth, vec = ' num2str(vec)])
    subplot(3,2,3); imagesc(RelErrorMap.SpiralDenoised_vs_SpirGroundTruth .* mask), colorbar, title('RelError SpiralDenoised vs SpiralGroundTruth')
    subplot(3,2,4); imagesc(RelErrorMap.SpiralSpice_vs_SpirGroundTruth .* mask), colorbar, title('RelError SpiralSpice vs SpiralGroundTruth')  
    subplot(3,2,5); imagesc(RelErrorMap.SpiralSpiceB0Corr_vs_SpirGroundTruth .* mask), colorbar, title('RelError SpiralSpiceB0Corr vs SpiralGroundTruth')  

    
    RelError
    
    figure; 
    subplot(3,2,1); imagesc(abs(SpSpice.GroundTruth.csi(:,:,:,vec)),[0 maxscale]), colorbar, title('Ground Truth')
    subplot(3,2,2); imagesc(abs(SpSpice.GroundTruth_Spi2Cart.csi(:,:,:,vec)),[0 maxscale]), colorbar, title(['GroundTruth Spi2Cart vec = ' num2str(vec)])
    subplot(3,2,3); imagesc(abs(SpSpice.Spi2Cart.csi(:,:,:,vec)),[0 maxscale]), colorbar, title('SpiralGriddedToCart') %
    subplot(3,2,4); imagesc(abs(SpSpice.Spi2Cart_Denoised.csi(:,:,:,vec)),[0 maxscale]), colorbar, title('SpiralGriddedToCart Denoised')
    subplot(3,2,5); imagesc(abs(SpSpice.Reco.csi(:,:,:,vec)),[0 maxscale]), colorbar, title('SpiralSpice')         
    subplot(3,2,6); imagesc(abs(SpSpice.Reco_B0Corr.csi(:,:,:,vec)),[0 maxscale]), colorbar, title('SpiralSpiceB0Corr')         
    

% end





%% Spectra


dt     = SpSpice.GroundTruth.Par.SpecDwellTimes(1)*1e-9;
BW     = 1/dt;
freq      = linspace(-BW/2,BW/2,SpSpice.GroundTruth.Par.MatSize(4));
ppm = 4.65-freq/SpSpice.GroundTruth.Par.LarmorFreq*10^6;

Voxels = {[27 26],[23 34],[30 44],[29 35],[40 39],[40 27]};
% Voxels = {[42 46],[27 26],[30 44],[30 36],[40 39],[40 27]};



figure;
SpirGroundTruth = SpSpice.GroundTruth_Spi2Cart.csi;
Reco_sFT = SpSpice.Spi2Cart.csi;
Reco_sFTDenoised = SpSpice.Spi2Cart_Denoised.csi;
Reco_SpiralSpice = SpSpice.Reco.csi;
Reco_SpiralSpiceB0Corr = SpSpice.Reco_B0Corr.csi;

for i = 1:numel(Voxels)
    subplot(3,2,i);
%     plot(ppm,real(fftshift(fft(squeeze(GroundTruth(Voxels{i}(1),Voxels{i}(2),1,:)))))/sqrt(SpSpice.GroundTruth.Par.MatSize(4)),'k'), 
    plot(ppm,real(fftshift(fft(squeeze(SpirGroundTruth(Voxels{i}(1),Voxels{i}(2),1,:)))))/sqrt(SpSpice.GroundTruth.Par.MatSize(4)),'k'),
    hold on
    plot(ppm,real(fftshift(fft(squeeze(Reco_sFT(Voxels{i}(1),Voxels{i}(2),1,:)))))/sqrt(SpSpice.GroundTruth.Par.MatSize(4)),'r'),
    plot(ppm,real(fftshift(fft(squeeze(Reco_sFTDenoised(Voxels{i}(1),Voxels{i}(2),1,:)))))/sqrt(SpSpice.GroundTruth.Par.MatSize(4)),'g'),
    plot(ppm,real(fftshift(fft(squeeze(Reco_SpiralSpice(Voxels{i}(1),Voxels{i}(2),1,:)))))/sqrt(SpSpice.GroundTruth.Par.MatSize(4)),'b')
    plot(ppm,real(fftshift(fft(squeeze(Reco_SpiralSpiceB0Corr(Voxels{i}(1),Voxels{i}(2),1,:)))))/sqrt(SpSpice.GroundTruth.Par.MatSize(4)),'m')
    hold off, legend('SpirGroundTruth','Reco sFT','Reco sFTDenoised', 'Reco SpiralSpice', 'Reco B0Spice')
    title(['InData x=' num2str(Voxels{i}(1)) ', y=' num2str(Voxels{i}(2))])
end


figure;
for i = 1:numel(Voxels)
    subplot(3,2,i);
    plot(ppm,abs(fftshift(fft(squeeze(SpirGroundTruth(Voxels{i}(1),Voxels{i}(2),1,:))))/sqrt(SpSpice.GroundTruth.Par.MatSize(4)) - fftshift(fft(squeeze(Reco_sFT(Voxels{i}(1),Voxels{i}(2),1,:))))/sqrt(SpSpice.GroundTruth.Par.MatSize(4))),'r')
    hold on
    plot(ppm,abs(fftshift(fft(squeeze(SpirGroundTruth(Voxels{i}(1),Voxels{i}(2),1,:))))/sqrt(SpSpice.GroundTruth.Par.MatSize(4)) - fftshift(fft(squeeze(Reco_sFTDenoised(Voxels{i}(1),Voxels{i}(2),1,:))))/sqrt(SpSpice.GroundTruth.Par.MatSize(4))),'g')
    plot(ppm,abs(fftshift(fft(squeeze(SpirGroundTruth(Voxels{i}(1),Voxels{i}(2),1,:))))/sqrt(SpSpice.GroundTruth.Par.MatSize(4)) - fftshift(fft(squeeze(Reco_SpiralSpice(Voxels{i}(1),Voxels{i}(2),1,:))))/sqrt(SpSpice.GroundTruth.Par.MatSize(4))),'b')
    plot(ppm,abs(fftshift(fft(squeeze(SpirGroundTruth(Voxels{i}(1),Voxels{i}(2),1,:))))/sqrt(SpSpice.GroundTruth.Par.MatSize(4)) - fftshift(fft(squeeze(Reco_SpiralSpiceB0Corr(Voxels{i}(1),Voxels{i}(2),1,:))))/sqrt(SpSpice.GroundTruth.Par.MatSize(4))),'m')
    hold off
    legend('Error Reco sFT','Error Reco sFTDenoised','Error Reco SpiralSpice','Error Reco B0Spice')
    title(['Error x=' num2str(Voxels{i}(1)) ', y=' num2str(Voxels{i}(2))])
end






%% 

Offset = 3;

figure; 
subplot(4,3,1+0*Offset), imagesc(abs(SpSpice.GroundTruth_Spi2Cart.csi(:,:,:,1))), colorbar, title('SpirGroundTruth Vec1')
subplot(4,3,1+1*Offset), imagesc(abs(SpSpice.Spi2Cart.csi(:,:,:,1))), colorbar, title('Spiral Vec1')
subplot(4,3,1+2*Offset), imagesc(abs(SpSpice.Reco.csi(:,:,:,1))), colorbar, title('SpiralSpice Vec1')
subplot(4,3,1+3*Offset), imagesc(abs(SpSpice.Reco_B0Corr.csi(:,:,:,1))), colorbar, title('SpiralSpiceB0Corr Vec1')

subplot(4,3,2+0*Offset), imagesc(abs(SpSpice.GroundTruth_Spi2Cart.csi(:,:,:,2))), colorbar, title('SpirGroundTruth Vec2')
subplot(4,3,2+1*Offset), imagesc(abs(SpSpice.Spi2Cart.csi(:,:,:,2))), colorbar, title('Spiral Vec2')
subplot(4,3,2+2*Offset), imagesc(abs(SpSpice.Reco.csi(:,:,:,2))), colorbar, title('SpiralSpice Vec2')
subplot(4,3,2+3*Offset), imagesc(abs(SpSpice.Reco_B0Corr.csi(:,:,:,2))), colorbar, title('SpiralSpiceB0Corr Vec2')

subplot(4,3,3+0*Offset), imagesc(abs(SpSpice.GroundTruth_Spi2Cart.csi(:,:,:,10))), colorbar, title('SpirGroundTruth Vec10')
subplot(4,3,3+1*Offset), imagesc(abs(SpSpice.Spi2Cart.csi(:,:,:,10))), colorbar, title('Spiral Vec10')
subplot(4,3,3+2*Offset), imagesc(abs(SpSpice.Reco.csi(:,:,:,10))), colorbar, title('SpiralSpice Vec10')
subplot(4,3,3+3*Offset), imagesc(abs(SpSpice.Reco_B0Corr.csi(:,:,:,10))), colorbar, title('SpiralSpiceB0Corr Vec10')


%% Metabolic Map Quantification

dt     = SpSpice.GroundTruth.Par.SpecDwellTimes(1)*1e-9;
BW     = 1/dt;
freq      = linspace(-BW/2,BW/2,SpSpice.GroundTruth.Par.MatSize(4));
ppm = 4.65-freq/SpSpice.GroundTruth.Par.LarmorFreq*10^6;
mask = SpSpice.GroundTruth.Masks.Vial{1};


CurData.GroundTruth = fftshift(fft(SpSpice.GroundTruth.csi,[],4),4)/sqrt(SpSpice.GroundTruth.Par.MatSize(4));
CurData.Reco_sFT = fftshift(fft(SpSpice.Spi2Cart.csi,[],4),4)/sqrt(SpSpice.GroundTruth.Par.MatSize(4));
CurData.Reco_sFTDenoised = fftshift(fft(SpSpice.Spi2Cart_Denoised.csi,[],4),4)/sqrt(SpSpice.GroundTruth.Par.MatSize(4));
CurData.Reco_SpiralSpice = fftshift(fft(SpSpice.Reco.csi,[],4),4)/sqrt(SpSpice.GroundTruth.Par.MatSize(4));
CurData.Reco_B0Spice = fftshift(fft(SpSpice.Reco_B0Corr.csi,[],4),4)/sqrt(SpSpice.GroundTruth.Par.MatSize(4));

SearchMetabos = {'Metabo', 'PeakPPM','IntegratePPMRange';'NAA', 2.01, 0.5; 'Cr', 3.02, 0.085; 'Cho',3.19, 0.085; 'H2O', 4.65, 0.5 };

RecoNames = transpose(fieldnames(CurData));
MetMapAbsError = cell([numel(RecoNames)+1 size(SearchMetabos,1)]);
MetMapAbsError(1,2:end) = {'NAA','Cr','Cho','H2O'};
MetMapAbsError(2:end,1) = RecoNames;

for CurRecoInd = 1:numel(RecoNames)
    CurReco = RecoNames{CurRecoInd};
    for MetNo = 2:size(SearchMetabos,1)
        [bla, PpmPointLow] = min(abs(ppm - (SearchMetabos{MetNo,2}+SearchMetabos{MetNo,3})));
        [bla, PpmPointHi] = min(abs(ppm - (SearchMetabos{MetNo,2}-SearchMetabos{MetNo,3})));
        MetMap.(CurReco).(SearchMetabos{MetNo,1}) = abs(sum(real(CurData.(CurReco)(:,:,:,PpmPointLow:PpmPointHi)),4));
        
        dummy = abs(MetMap.(CurReco).(SearchMetabos{MetNo,1}) - MetMap.GroundTruth.(SearchMetabos{MetNo,1}));
%         dummy = dummy ./ abs(MetMap.GroundTruth.(SearchMetabos{MetNo,1}));
        dummy(mask == 0) = NaN;
        MetMapAbsError{CurRecoInd+1,MetNo} = nanmean_own(dummy(:));
    end
end
clear bla CurData PpmPointLow PpmPointHi dummy


SubPlotRows = size(SearchMetabos,1)-1;
SubPlotCols = numel(fields(MetMap));
figure; 
% for CurRecoInd = 1:numel(RecoNames)
%     CurReco = RecoNames{CurRecoInd};
%     for MetNo = 1:size(SearchMetabos,1)-1
%         subplot(SubPlotRows,SubPlotCols,(CurRecoInd-1)*SubPlotCols + MetNo)
%         imagesc(MetMap.(CurReco).(SearchMetabos{MetNo+1,1}))
%         title([CurReco ', ' SearchMetabos{MetNo+1,1}], 'Interpreter','none')
%     end
% end
[ha, pos] = tight_subplot(SubPlotCols,SubPlotRows,0.05,0.05,0);

ii = 0;
for CurRecoInd = 1:SubPlotCols
    CurReco = RecoNames{CurRecoInd};
    for MetNo = 1:SubPlotRows
        ii = ii+1;
        axes(ha(ii)); imagesc(MetMap.(CurReco).(SearchMetabos{MetNo+1,1}))
        title([CurReco ', ' SearchMetabos{MetNo+1,1}], 'Interpreter','none')
    end
end



MetMapAbsError
clearvars -except MetMapAbsError MetMap SpSpice RelErrorMap *_flag


