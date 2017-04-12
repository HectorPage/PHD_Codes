function SO_search(param_one, tau, param_two)

%This program performs a 3D parameter search for phiE1I1 and phiI1E1 in the
%amputated oscillator for a specific RC range.
%Expects parameter ranges given as arrays.
%Starts in parent folder 'parentpath', in which you need to put the .exe 'gaussian_oscillator'.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hector JI Page Date: 30/11/11         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for adx = 1:numel(param_two)
    dirstringT = [num2str(param_two(adx))];
    parentpath = ['~/video_conflict/ff_plasticity/_moving_rat/full_COMB_models/toy_model/MultiThread/_tanh/_RC_effect/nonrandom'];
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
            tier_3_path = [tier_2_path, '/phitest', dirstringI];
		
            if ~exist(tier_3_path, 'dir')
            
            mkdir(tier_3_path);
            copyfile('~/video_conflict/ff_plasticity/_moving_rat/full_COMB_models/toy_model/MultiThread/_tanh/_RC_effect/nonrandom/RC_effect', tier_3_path);
            cd(tier_3_path);
            
        
            fid = fopen('RCparams.dat', 'w+');  
            fprintf(fid, 'time = 300.5\n');
            fprintf(fid, 'training time = 298.5\n');
            fprintf(fid, 'E1 cells = 500\n');
            fprintf(fid, 'RC1 connections = 500\n');    
            fprintf(fid, 'phi RC1 = 60.0\n');
            fprintf(fid, ['phi RC1 test = ', dirstringI, '\n']);
            fprintf(fid, 'tau E1 = 0.001\n');
            fprintf(fid, 'alpha E1 = 0.0\n');
            fprintf(fid, 'beta E1 = 0.3\n');
            fprintf(fid, 'timestep size = 0.0001\n');
            fprintf(fid, 'external inhibition = 50.0\n');
            fprintf(fid, 'global inhibition = 0.01\n');
            fprintf(fid, 'input strength = 70.0\n'); 
            fprintf(fid, 'input location = 180.0\n');
            fprintf(fid, 'input sigma = 30.0\n');
            fprintf(fid, 'learning rate = 0.01\n');
            fprintf(fid, 'normalise = 0\n');
            fprintf(fid, 'velocity = 180.0\n');
            fprintf(fid, 'conduction delay = 0.01\n');
            fprintf(fid, 'parallel threads = 10\n');
            fprintf(fid, 'sigmoid = 0\n');
            fprintf(fid, 'symmetrical testing = 0\n');
            fprintf(fid, 'symmetrical sigma = 20.0\n');
            fprintf(fid, 'symmetrical strength = 0\n');
            fprintf(fid, 'start length = 0\n');
            fprintf(fid, 'random weights = 0');
           
            
		
            fclose(fid);
        
            system(['./RC_effect -fname RCparams.dat']);
            
            plot_toy(300.5, 0.0001, 500, 'save', 0);
            toy_weights_plots([50:50:500], 10, 1, 1, 500, 180.0, 0.01);
            offset_evo(500, 300.5, 0.0001, 0.00001, 1, 0);
            
            
            delete('RC_effect');
           
            
            
            cd(tier_2_path);
            
            else
                cd(tier_3_path);
                plot_toy(300.5, 0.0001, 500, 'save', 0);
                toy_weights_plots([50:50:500], 10, 1, 1, 500, 180.0, 0.01);          
                offset_evo(500, 300.5, 0.0001, 0.00001, 1, 0);
            
            end
            
           
            
            cd(tier_2_path);
            
        end
    
        cd(tier_1_path); 
   
    end
    
    cd(parentpath);

end

end
