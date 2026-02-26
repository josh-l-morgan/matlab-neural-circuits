function [all_results, max_workers]=bundle_and_distribute_viafile(job_function, distributed_args, max_workers, job_manager, path_dependencies, compile_all_results, varargin)

EMroot='/home/viren/EM';

permutation=randperm(length(distributed_args));

% split job into bundles
bundle_size=ceil(length(distributed_args)/max_workers);
distributed_bundle=cell(max_workers, 1);

arg_num=0;
for i=1:max_workers

	bundle=[];
	
	for j=1:bundle_size
		arg_num=arg_num+1;
		if(arg_num<=length(distributed_args))
			
			arguments=distributed_args{permutation(arg_num)};
			arguments{end+1}=job_function;			
			bundle{j}=arguments;
		end
	end
	
	distributed_bundle{i}={bundle};
	
end

if(nargin>6)
	job_data=varargin{1};

	rand_n=num2str(ceil(rand*100000));
	jobdata_fname=get_temp_filename(['/home/viren/process_networks/bundle_temp/bundle_temp_jobdata_', rand_n, '_'], '.mat');
	save(jobdata_fname, 'job_data');
end

display('writing out distributed arguments to file...');
for i=1:length(distributed_bundle)
	rand_n=num2str(ceil(rand*100000));
	bundle_fname=get_temp_filename(['/home/viren/process_networks/bundle_temp/bundle_temp_', rand_n, '_'], '.mat');
	bundle=distributed_bundle{i};
	save(bundle_fname, 'bundle');
	if(nargin>6)
		distributed_bundle{i}={bundle_fname, jobdata_fname};
	else
		distributed_bundle{i}={bundle_fname};
	end
end


display('submitting job to DCE...');
sched=get_sched(job_manager);
j=sched.createJob();
evaldiff=createTask(j, @bundle_and_distribute_worker_viafile, 0, distributed_bundle);
set(j, 'MaximumNumberOfWorkers', max_workers);
set(j, 'PathDependencies',path_dependencies);
set(evaldiff,'CaptureCommandWindowOutput',true);
jobID=j.ID;
display(['distributed processing on job id ', num2str(jobID),'.']);
% run the job
submit(j);
tic, waitForState(j);, toc
errmsgs=get(evaldiff,{'ErrorMessage'});
msgs=get(evaldiff,{'CommandWindowOutput'});

% combine all bundled results
results=get(evaldiff,'OutputArguments');

if(compile_all_results)
	
	all_results=cell(length(distributed_args), 1);
	arg_num=0;
	
	try
	
	for i=1:max_workers
		bundle_fname=[distributed_bundle{i}{1}, '_return.mat'];
		bundle_results=load(bundle_fname);
		bundle_results=bundle_results.bundle_results;
		for k=1:bundle_size
			arg_num=arg_num+1;
			if(arg_num<=length(distributed_args))
				all_results{permutation(arg_num)}=bundle_results{k};			
			end
		end	
	end
	
	catch
		keyboard
	end
else

	all_results=distributed_bundle;
end

j.destroy;