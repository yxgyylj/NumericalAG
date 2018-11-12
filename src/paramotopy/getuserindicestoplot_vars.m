% [nplot] = getuserindicestoplot_vars(info)
%
% gets input from the user on how many variables to plot, and which ones.
%
% input: info, a struct generated by parse_input()
%        info has the following data members:
%              filename: 'alicia'
%                numeqn: 9
%           numvargroup: 1
%              numparam: 15
%              numconst: 0
%                   eqn: {1x9 cell}
%                  vars: {1x9 cell}
%           userdefined: 1
%            paramnames: {15x1 cell}
%     userparameterfile: 'spacesample_15_10000'
%         supplementary: {[1x69 char]}
%                numvar: 9
%              mc_lines: 10000
%
% output: nplot, a struct with indices and a variable number.
%            nplot has the following fields:
%       userparams: 3
%          indices: [integer]
%
% daniel brake
% colorado state university
% mathematics
% 2013
% danielthebrake@gmail.com

function [nplot] = getuserindicestoplot_vars(info)

display(info.vars)
nplot.uservars = input('there are more than 3 variables, and simple no-input plotting supported\nonly for three or less...\ntherefore, how many variables (<=3) would you like to plot?\nnote:can use 4, if 4th is color\n');
for ii = 1:info.numvar
	display(sprintf([num2str(ii) ': ' info.vars{ii}]));
end

nplot.indices = zeros(nplot.uservars,1);
for ii = 1:nplot.uservars
	nplot.indices(ii) = input(['which variable in the ' num2str(ii) 'th position?']);
end

end

