function [all_slopes] =  get_tf_slopes(x_axis,A,varargin)

if ~isempty(varargin)
    sig_threshold = varargin{1};
    for i_c = 1:size(A,2)
        sig_mask(:,:,i_c) = A{1,i_c} > repmat(sig_threshold{1,i_c}(:,1),1,size(A{1,i_c},2));
    end
end



all_slopes = nan([size(A{1,1})]);
for i_freq = 1:size(A{1,1},1)
    for i_t = 1:size(A{1,1},2)
        curr_pol = polyfit(x_axis,cellfun(@(x) x(i_freq,i_t),A),1);
        all_slopes(i_freq,i_t) = curr_pol(1);
    end
end

% all_change_pts= nan([size(A{1,1})]);
% all_residuals = nan([size(A{1,1})]);
% for i_freq = 1:size(A{1,1},1)
%     for i_t = 1:size(A{1,1},2)
%     [ipt residual] = findchangepts(cellfun(@(x) x(i_freq,i_t),A));
%     all_residuals(i_freq,i_t) = residual;
%     all_change_pts(i_freq,i_t) =ipt;
%     end
% end