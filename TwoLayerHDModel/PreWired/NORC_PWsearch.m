function NORC_PWsearch(param_one, tau, param_two, SO)

%This program performs a 3D parameter search for phiE1I1 and phiI1E1 in the
%amputated oscillator for a specific RC range.
%Expects parameter ranges given as arrays.
%Starts in parent folder 'parentpath', in which you need to put the .exe 'gaussian_oscillator'.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hector JI Page Date: 30/11/11         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for adx = 1:numel(param_two)
    dirstringT = [num2str(param_two(adx))];
    parentpath = ['~/video_conflict/ff_plasticity/_moving_rat/full_COMB_models/no_RC/single_layer/MultiThread/BothDelays/parameter_search/1/1/_new_tau0.001'];
    tier_1_path = [parentpath,'/',dirstringT];
    
    if ~exist(tier_1_path, 'dir')
        mkdir(tier_1_path);
    end
  
    cd(tier_1_path);       

    for idx = 1:numel(param_one)
        dirstringO = [num2str(param_one(idx))];
        tier_2_path = [tier_1_path,'/',dirstringO];
    
        if exist(tier_2_path, 'dir')        
            cd(tier_2_path); 
        else
            mkdir(tier_2_path);
            cd(tier_2_path);
        
        end
    
        for jdx = 1:numel(tau)
    
            dirstringI = [num2str(tau(jdx))];
            tier_3_path = [tier_2_path, '/delay', dirstringI];
		
            if ~exist(tier_3_path, 'dir')
            
            mkdir(tier_3_path);
            copyfile('~/video_conflict/ff_plasticity/_moving_rat/full_COMB_models/no_RC/single_layer/MultiThread/BothDelays/norcV2', tier_3_path);
            cd(tier_3_path);
            
        
            fid = fopen('norc_parameters.dat', 'w+');  
            fprintf(fid, 'time = 4.1\n');
            fprintf(fid, 'pause time = 1.0\n');
            fprintf(fid, 'rotation time = 2.0\n');
            fprintf(fid, 'input time = 0.1\n');
            fprintf(fid, 'HD cells = 500\n');
            fprintf(fid, 'COMB cells = 1000\n');
            fprintf(fid, 'ROT cells = 1\n');
            fprintf(fid, 'NOROT cells = 1\n'); 
            fprintf(fid, 'phi HDCOMB = 700.0\n');
            fprintf(fid, 'phi COMBHD = 4500.0\n');   
            fprintf(fid, 'phi ROT = 80.0\n');
            fprintf(fid, 'phi NOROT = 80.0\n');
            fprintf(fid, 'tau HD = 0.001\n');  
            fprintf(fid, 'tau COMB = 0.001\n'); 
            fprintf(fid, 'alpha HD = 0.0\n');
            fprintf(fid, 'beta HD = 0.2\n');
            fprintf(fid, 'alpha COMB = 16.0\n');
            fprintf(fid, 'beta COMB = 0.3\n');
            fprintf(fid, 'timestep size = 0.00001\n');        %Originally 0.00001
            fprintf(fid, 'external inhibition = 0.0\n');        %originally 0.0
            fprintf(fid, 'global inhibition = 0.2\n');          %originally 0.2
            fprintf(fid, 'COMB inhibition = 0.35\n');           %originally 0.35
            fprintf(fid, 'initialiser strength = 2.0\n');       %originally 2.0
            fprintf(fid, 'initialiser location = 180.0\n');
            fprintf(fid, 'initialiser sigma = 20.0\n');
            fprintf(fid, 'HD to COMB weight width = 20.0\n');   
            fprintf(fid, 'COMB to HD weight width = 20.0\n');   
            fprintf(fid, ['conduction delay = ',dirstringI ,'\n']);
            fprintf(fid, 'velocity = 180.0\n');
            fprintf(fid, 'threads = 10');
           
            
		
            fclose(fid);
        
            system(['./norcV2 -fname norc_parameters.dat']);
            
            if(SO)
            norc_plot(4.1, 0.00001, 500, 'save');
            norc_oscil(4.1, 0.00001, 500, 'save');
            %norc_weights([250:500:1000], 2, 1, 1, 500, 1000, 180.0, 0.05, 1);
            else
            norc_plot(4.1, 0.00001, 500, 'save');
            norc_oscil(4.1, 0.00001, 500, 'save');
            %norc_weights([250:500:1000], 2, 1, 1, 500, 1000, 180.0, 0.05, 0);
            end
            
            delete('norcV2');
            delete('COMBActivations.bdat');
            delete('HDActivations.bdat');
            delete('_COMBrates.bdat');
            delete('_HDrates.bdat');
            delete('_speeds.bdat');
            delete('HDrates.bdat');
            delete('COMBrates.bdat');
            delete('COMBtoHDweights.bdat');
            delete('HDtoCOMBweights.bdat');
            delete('Input_Location1.bdat');
            delete('NOROTRates.bdat');
            delete('ROTRates.bdat');
            delete('Pvector.dat.bdat');
            
            
            cd(tier_2_path);
            
            else
                cd(tier_3_path);
                 
            if(SO)
            norc_plot(4.1, 0.00001, 500, 'save');
            norc_oscil(4.1, 0.00001, 500, 'save');
            %norc_weights([250:500:1000], 2, 1, 1, 500, 1000, 180.0, 0.05, 1);
            else
            norc_plot(4.1, 0.00001, 500, 'save');
            norc_oscil(4.1, 0.00001, 500, 'save');
           %norc_weights([250:500:1000], 2, 1, 1, 500, 1000, 180.0, 0.05, 0);
            end
            end
            
           
            
            cd(tier_2_path);
            
        end
    
        cd(tier_1_path); 
   
    end
    
    cd(parentpath);

end

end
