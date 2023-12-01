% This code contains the main body of calculations for homogenization of 
% plastic REV

clc
clear
format long g
tic;

% load the input properties
[LinkMicro,LoadParam,gparam,ParamMicro,PropMicro]=InputFile();

% initialize different variables
[StrainMic,StressMic,StrainMac,StressMac,ParamMicro]= ...
    Initialize(ParamMicro,gparam);

% caculate localization tensors and the homogenized elasticity
% as we neglect change in configurational properties, this is performed
% once before the start of the analysis
[TensorsMicro,ParamMicro,CmatLocalMicro,PropMicro]= ...
    Concentration(ParamMicro,PropMicro,LinkMicro,gparam);

% update output variable
[out(1).strain]=Record(StrainMac);
[out(1).stress]=Record(StressMac);

for istage=1:size(LoadParam,2)

    FacStep=0; % for tracing the total number of strain increments applied to the sample
    for istep=1:LoadParam(istage).nstep

        if (istep==1)
            formatspec='-----------------stage number-----------------= %d\n';
            fprintf(formatspec,istage)
        end

        if (istep==19)
            kashk=1;
        end
        
        formatspec='--------------step number--------------= %d\n';
        fprintf(formatspec,istep)
        pause(0.5)

        UnloadIndic=0; % unloading step indicator (activated after a desired step)
        if (istep > 50e10)
            UnloadIndic=1;
        end
        
        DelFacStep=1;
        if (UnloadIndic==1)
            DelFacStep=-1;
        end
        FacStep=FacStep+DelFacStep;

        % total current macro strain in the sample (prescribed)
        StrainStage=Convert(FacStep*LoadParam(istage).StrainStep,1);
        StrainMac.cur=StrainStage+StrainMac.stage;
        StressStage=Convert(FacStep*LoadParam(istage).StressStep,1);
        StressMac.cur=StressStage+StressMac.stage;

        % calculate the full macro strain tensor considering prescribed
        % stresses
        [StrainMac]=SolveMacroStrain(ParamMicro,StrainMac,StressMac,LoadParam(istage));
        
        % to localize strain, do return mapping which includes updating
        % micro strains and stresses
        [YieldMicro,ParamMicro,StrainMic,ReturnIter,StrainMac,StressMic,StressMac]= ...
            ReturnMap(ParamMicro,StrainMic,LinkMicro,TensorsMicro, ...
            StrainMac,CmatLocalMicro,PropMicro,gparam,StressMac, ...
            LoadParam(istage),StressMic);
         
        % add cross-checks here

        % update output variable
        [out(istep+1).strain]=Record(StrainMac);
        [out(istep+1).stress]=Record(StressMac);
        
        % initialize next step
        StrainMac.plpre=StrainMac.plcur;
        StrainMac.pre=StrainMac.cur;
        StressMac.pre=StressMac.cur;
        for ipoin=1:ParamMicro.npoin
            StrainMic(ipoin).plpre=StrainMic(ipoin).plcur;
            StrainMic(ipoin).pre=StrainMic(ipoin).cur;
            StressMic(ipoin).pre=StressMic(ipoin).cur;
        end
        
    end % step

    % initialize next stage
    StrainMac.stage=StrainMac.cur;
    StressMac.stage=StrainMac.cur;

end % stage

save ('AllVars.mat')