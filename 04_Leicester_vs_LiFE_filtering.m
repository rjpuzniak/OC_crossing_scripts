%% Initialization

% Add required folders to path
addpath(genpath('/home/auguser2016/Software/spm12'))
addpath(genpath('/home/auguser2016/Software/vistasoft'))
addpath(genpath('/home/auguser2016/Software/mba'))
addpath(genpath('/home/auguser2016/Software/encode'))
addpath(genpath('/home/auguser2016/Software/AFQ'))

% Prepare list of subjects
cd /home/auguser2016/Projects/0001_Chiasm/
subjects=dir

% analyzed variant
all_variants={'CSD_seed';'CSD_selected';'DT_seed';'DT_selected'}%;'DT2_seed';'DT2_selected'}
%all_variants={'DT2_selected'}

% Leicester or LiFE
all_filterings={false, true}
%all_filterings={false}

for xx=1:length(all_variants)
    for yy=1:length(all_filterings)
        
        variant=all_variants{xx};
        filtering=all_filterings{yy};

        % Initialize tables storing results
        no_er = array2table(nan(2,size(subjects,1)-4));
        no_er.Properties.VariableNames = {subjects(5:end).name};
        no_er.Properties.RowNames = {'Fascicles Remaining',' Mean Error'};
        no_er

        con_ips = array2table(nan(4,size(subjects,1)-4));
        con_ips.Properties.VariableNames = {subjects(5:end).name};
        con_ips.Properties.RowNames = {'Right to Right',' Right to Left','Left to Right','Left to Left'};
        con_ips

        sanity_check = array2table(nan(1,size(subjects,1)-4))
        sanity_check.Properties.VariableNames = {subjects(5:end).name};
        sanity_check

        cd('/home/auguser2016/Projects/Leicester_vs_LiFE')
        current = pwd;

        %% Iterate through all subjects
        %iterate=[5, 7:size(subjects,1)]
        for i=5:size(subjects,1)%i=5 failed for ce04

            subj=subjects(i).name;

            % Define relevant folders
            auxilary_data='/home/auguser2016/dMRI_DATA/PREPROCESSED_DATA/';
            t1File        = fullfile(strcat(auxilary_data,'Anatomies_ACPC_Aligned/',subj,'_t1_acpc.nii.gz'));

            specific_data=strcat('/home/auguser2016/Projects/Leicester_vs_LiFE/Leicester_',variant,'/',subj,'/');

            dwiFile       = fullfile(strcat('/home/auguser2016/Projects/0001_Chiasm/',subj,'/',subj,'_150mm_clean_aligned_trilin_noMEC.nii.gz'));
            dwiFileRepeat = fullfile(strcat('/home/auguser2016/Projects/0001_Chiasm/',subj,'/',subj,'_150mm_clean_aligned_trilin_noMEC.nii.gz'));

            % Build the file names for the diffusion data, the anatomical MRI.
            fgFileName = fullfile(strcat(specific_data,subj,'_ACT.tck'))

            % The final connectome and data astructure will be saved with this name:
            feFileName    = strcat(subj,'_postlife_chiasm');

            % Discretization parameters
            L = 360; 
            Niter = 500;

            % check wwhether everything is all right           
            [status, result]=system(strcat('tckstats',32,fgFileName,32,'-output count -quiet'))
            if (result==0) 
                continue
            end
             
            % initial - fe, fg and track_name will be overwritten in case of
            % filtering
            
            try
                fe = feConnectomeInit(dwiFile,fgFileName,feFileName,[],dwiFileRepeat,t1File,L,[0,1]);
            catch
                continue % but do nothing with error, just skip this iteration of loop
            end
            
            fg = feGet(fe,'fibers acpc');
            track_name='noLiFE'

            if istrue(filtering)

                % Initializing structure and fitting model
                fe = feSet(fe,'fit',feFitModel(feGet(fe,'model'),feGet(fe,'dsigdemeaned'),'bbnnls',Niter,'preconditioner'));

                % Extract non zero weighted fascicles
                fg = feGet(fe,'fibers acpc');
                w = feGet(fe,'fiber weights');
                positive_w = w > 0;
                
                
                try fg = fgExtract(fg, positive_w, 'keep');
                catch
                    continue
                end
                    
                % Get mean error and write results to table
                mean_error = mean(feGet(fe,'vox rmse'))
                no_er.(subj) = [sum(positive_w); mean_error]

                track_name='postLiFE'

            end

            % Save results (fibers, images) to separate folder
            cd (specific_data);  
            fgWrite(fg,strcat(subj,'_', track_name, '_chiasm.tck'),'tck')

            % define a folder for filtered and non-filtered results
            if istrue(filtering)
                output_results_folder=strcat(specific_data,subj,'_','postLiFE')
            else
                output_results_folder=strcat(specific_data,subj,'_','noLiFE')
            end
            mkdir(output_results_folder)

            % write segmentation results to given folder
            system(sprintf('source /home/auguser2016/Projects/Leicester_vs_LiFE/07_connecting_rois.sh %s %s %s', subj, strcat(pwd,'/',subj,'_', track_name, '_chiasm.tck'), output_results_folder))  

            % Check if the requested file is there. If not write 0 in required place

            [~,l2l]=system(strcat('tckstats',32,output_results_folder,'/',subj,'_left_to_left.tck -output count -quiet'),'-echo')
            [~,l2r]=system(strcat('tckstats',32,output_results_folder,'/',subj,'_left_to_right.tck -output count -quiet'),'-echo')
            [~,r2l]=system(strcat('tckstats',32,output_results_folder,'/',subj,'_right_to_left.tck -output count -quiet'),'-echo')
            [~,r2r]=system(strcat('tckstats',32,output_results_folder,'/',subj,'_right_to_right.tck -output count -quiet'),'-echo')

            if(size(l2l,2)==87)
                l2l='0'
            end

            if(size(l2r,2)==87)
                l2r='0'
            end

            if(size(r2l,2)==87)
                r2l='0'
            end

            if(size(r2r,2)==87)
                r2r='0'
            end

            % Sanity check
            sanity_check.(subj) = -size(fg.fibers,1) + str2num(l2l) + str2num(l2r) + str2num(r2l) + str2num(r2r)

            % Write results to our table
            %con_ips.(subj)=[str2num(r2r);str2num(r2l);str2num(l2r);str2num(l2l)]

            con_ips.(subj)('Right to Right')=str2num(r2r)
            con_ips.(subj)('Right to Left')=str2num(r2l)
            con_ips.(subj)('Left to Right')=str2num(l2r)
            con_ips.(subj)('Left to Left')=str2num(l2l)

        end 

        %finally save the table with the right name to relevant location
        save(strcat('/home/auguser2016/Projects/Leicester_vs_LiFE/', variant,'_',track_name),'con_ips')

    end
    
end