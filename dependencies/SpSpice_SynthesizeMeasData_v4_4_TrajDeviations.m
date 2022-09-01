function Output = SpSpice_SynthesizeMeasData_v4_4_TrajDeviations(Transj,Input,V,SamplingOperator,B0Shiftmap,DCFPreG,sft2_Oper,TiltTrajMat)
%
% Spice_SynthesizeMeasData 
%
% This function was written by Bernhard Strasser, April 2018.
%
%
% 
%
%
% Output = Spice_SynthesizeMeasData(U,V,B0Shiftmap,SamplingOperator)
%
% Input: 
% -         U                    ...     
% Output:
% -         Output                         ...     



% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!

% Further remarks: 




%% 0. Preparations

if(~exist('Transj','var') || isempty(Transj))
   Transj = 'NoTransj'; 
end



%%


if(strcmpi(Transj,'NoTransj') || strcmpi(Transj,'notransp')) % Going from UTS to NonCart k-space data

    
    Input = reshape(Input,[size(B0Shiftmap,1) size(V,1)]);

    % Synthesize MRSI Data
    Output = Input*V;
    
    % Apply B0-Effect
    Output = Output .* B0Shiftmap;

    % Inverse NUFFT (i-Space --> NonCart k-Space)
    Output = reshape(Output,[size(sft2_Oper,2) numel(Output)/size(sft2_Oper,2)]);
    Output = sft2_Oper * Output;
    
%     dummy = Output(:,1);
    Output = fft(fftshift(TiltTrajMat.*fftshift(ifft(Output,[],2),2),2),[],2);
%     Output(:,1) = dummy;

    
    Output = SamplingOperator .* Output;
    
    
    
    

    
    

else    % Going from NonCart k-space data to UTS

    
    Output = reshape(Input,size(SamplingOperator));
    Output = SamplingOperator .* Output;
    
%     dummy = Output(:,1);
    Output = fft(fftshift(conj(TiltTrajMat).*fftshift(ifft(Output,[],2),2),2),[],2);
%     Output(:,1) = dummy;    
    
    % NUFFT (NonCart k-Space --> i-Space)
    Output = Output .* myrepmat(DCFPreG,size(Output)); % Apply PrGridding Density Comp
	Output = sft2_Oper' * Output * size(sft2_Oper,2);
    
    
    % Undo B0-Effect
    Output = reshape(Output,size(B0Shiftmap));
    Output = Output .* conj(B0Shiftmap);

    % Synthesize MRSI Data
    Output = Output * V';

    Output = reshape(Output,[numel(Output) 1]);
    
    
    
    
    

    
    
    
    
    
end


