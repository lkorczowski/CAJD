function options = gp_getOpt(options,opt_dflt)

% deal with an input options structure and fill it according to a template
% opt_dflt.
%
% Inputs:
%   - options:  options structure input
%   - opt_dflt: template, default options structure of the function
%
% Output:
%   - options: final options structure used by the function
%
% *** History: 05-Jan-2015
% *** Author: Louis KORCZOWSKI, Florent BOUCHARD, GIPSA-Lab, 2015


% read the acceptable names
optionNames = fieldnames(opt_dflt);
% deal whith input
if isempty(options)
    warning('default options')
else
    nameArgs=fieldnames(options);
    nArgs = length(nameArgs);
    
    for indO = 1: nArgs
        inpName = nameArgs{indO};
        if any(strcmp(inpName,optionNames))
            opt_dflt.(inpName) = options.(inpName);
        else
            warning('%s is not a recognized parameter name',inpName)
        end
    end
end
% fill output structure
options = opt_dflt;
end