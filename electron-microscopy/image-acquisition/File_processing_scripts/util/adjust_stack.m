function [stack_adjusted]=adjust_stack(stack, varargin)

if(nargin>1)
	bb=varargin{1};
	stack_norm=stack(bb(1,1):bb(1,2), bb(2,1):bb(2,2), bb(3,1):bb(3,2));
else
	stack_norm=stack;
end

minstack=min(stack_norm(:));
if(minstack<0)
	stack=stack+abs(minstack);
end

stack=stack./max(stack_norm(:));
stack_adjusted=stack;

for i=1:size(stack,3)
	stack_adjusted(:,:,i)=imadjust(stack(:,:,i));
end