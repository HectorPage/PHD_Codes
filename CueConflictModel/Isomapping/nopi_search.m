function nopi_search(param_one, location, param_two)

%This program performs a 3D parameter search for phiE1I1 and phiI1E1 in the
%amputated oscillator for a specific RC range.
%Expects parameter ranges given as arrays.
%Starts in parent folder 'parentpath', in which you need to put the .exe 'gaussian_oscillator'.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Hector JI Page Date: 30/11/11         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for adx = 1:numel(param_two)
    dirstringT = [num2str(param_two(adx))];
    parentpath = ['~/video_conflict/ff_plasticity/_moving_rat/_no_PI/parameter_search/MultiThread'];
    tier_1_path = [parentpath,'/deg',dirstringT];
    
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
    
        for jdx = 1:numel(location)
    
            dirstringI = [num2str(location(jdx))];
            tier_3_path = [tier_2_path, '/', dirstringI];
		
            if ~exist(tier_3_path, 'dir')
            
            mkdir(tier_3_path);
            copyfile('~/video_conflict/ff_plasticity/_moving_rat/_no_PI/parameter_search/MultiThread/nopi', tier_3_path);
            cd(tier_3_path);
        
            fid = fopen('nopi_parameters.dat', 'w+');  
            fprintf(fid, 'time = 25.0\n');
            fprintf(fid, 'conflict time = 2.0\n');
            fprintf(fid, 'input time = 0.1\n');
            fprintf(fid, 'integration time = 22.0\n');
            fprintf(fid, 'HD cells = 500\n');
            fprintf(fid, 'RC1 connections = 500\n');    
            fprintf(fid, 'phi RC1 = 200.0\n');		
            fprintf(fid, 'tau HD = 0.001\n');   
            fprintf(fid, 'alpha HD = 0.0\n');
            fprintf(fid, 'beta HD = 0.3\n');
            fprintf(fid, 'timestep size = 0.0001\n');
            fprintf(fid, 'external inhibition = 0.0\n');
            fprintf(fid, 'global inhibition = 0.3\n');
            fprintf(fid, 'pi strength = 180\n'); 
            fprintf(fid, 'pi location = 180.0\n');
            fprintf(fid, 'pi sigma = 20.0\n');   
            fprintf(fid, 'RC weight width = 35.0\n');
            fprintf(fid, ['visring weight width = ', dirstringO, '\n']);
            fprintf(fid, 'phi visring = 6.5\n');
            fprintf(fid, 'visring strength = 2.0\n');
            fprintf(fid, ['visring location = ',dirstringI, '\n']);
            fprintf(fid, 'visring sigma = 50.0\n');
            fprintf(fid, 'learning rate = 0.05\n');
            fprintf(fid, 'normalise = 1\n');
            fprintf(fid, ['velocity = ',dirstringT,'\n']);
            fprintf(fid, 'parallel threads = 10');
            
            
		
            fclose(fid);
        
            system(['./nopi -fname nopi_parameters.dat']);
            
            plot_ff_touretzky(25.0, 0.0001, 500, 'save', 'Y', 'hd');
            
            delete('nopi');
           
            
            
            cd(tier_2_path);
        
            end
            
            cd(tier_2_path);
            
        end
    
        cd(tier_1_path); 
   
    end
    
    cd(parentpath);

end

end
