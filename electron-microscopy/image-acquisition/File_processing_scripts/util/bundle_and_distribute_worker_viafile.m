function []=bundle_and_distribute_worker(bundle_fname, varargin)

addpath /home/viren/EM/util/
addpath /home/viren/EM/orientation/
addpath /home/viren/EM/nn/hdf5/
addpath /home/viren/EM/analysis/
addpath /home/viren/EM/segmentation/
addpath /home/viren/EM/supergraph/iterated/
addpath /home/viren/EM/supergraph/
addpath /home/viren/EM/splits/watershed/tree/
addpath /home/viren/EM/splits/watershed/

bundle=load(bundle_fname, 'bundle');
bundle=bundle.bundle{1};
bundle_results=[];

if(nargin > 1)
	job_data=load(varargin{1});
	job_data=job_data.job_data;
	
end
for i=1:length(bundle)
	args=bundle{i};
	job_f=args{end};
	bundle_results{i}=cell(1,nargout(job_f));
	if(nargin > 1)
		[bundle_results{i}{:}]=job_f(args{1:end-1}, job_data{:});
	else
		[bundle_results{i}{:}]=job_f(args{1:end-1});
	end
end
 
rand_n=ceil(rand*10000);
bundle_fname=[bundle_fname, '_return.mat'];
save(bundle_fname, 'bundle_results');
