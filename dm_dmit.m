function [eval, EVect, EVect2, Components, S, Ddiff] = dm_dmit(X, varargin)
[S, D] = dm_simmat(X, varargin{:});

if nargout >= 6 % we need to find the diffusion distance in EVect2 space
    [eval, EVect, EVect2, Components, Ddiff] = dm_dmsim(S, 'Dist', D, varargin{:});
else
    [eval, EVect, EVect2, Components] = dm_dmsim(S, 'Dist', D, varargin{:});
end



end
