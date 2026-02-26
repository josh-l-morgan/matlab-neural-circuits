function [bundle_results]=bundle_and_distribute_worker(bundle)

javaaddpath /home/viren/EM/splits/watershed/java
addpath /home/viren/EM/util/
addpath /home/viren/EM/orientation/
addpath /home/viren/EM/nn/hdf5/
addpath /home/viren/EM/analysis/
addpath /home/viren/EM/segmentation/
addpath /home/viren/EM/supergraph/from_skeletons/
addpath /home/viren/EM/supergraph/
addpath /home/viren/EM/splits/watershed/tree/
addpath /home/viren/EM/splits/watershed/

bundle_results=[];
for i=1:length(bundle)
	args=bundle{i};
	job_f=args{end};
	bundle_results{i}=cell(1,nargout(job_f));
	[bundle_results{i}{:}]=job_f(args{1:end-1});
end
