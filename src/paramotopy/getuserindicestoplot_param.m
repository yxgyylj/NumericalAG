% [nplot] = getuserindicestoplot_param(info)
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



function [nplot] = getuserindicestoplot_param(info)

display(info.paramnames)

nplot.userparams = -1;

while or(nplot.userparams>3,nplot.userparams<1)
	tmp = input('there are more than 3 parameters, and simple no-input plotting supported\nonly for three or less...\ntherefore, how many parameters (<=3) would you like to plot?\n');
	while isempty(tmp)
		tmp = input('');
	end
	nplot.userparams = tmp;
end


for ii = 1:info.numparam
	if ii<10
		tmp = ' ';
	else
		tmp = '';
	end
	display(sprintf([tmp num2str(ii) ': ' info.paramnames{ii}]));
end

nplot.indices = zeros(nplot.userparams,1);
for ii = 1:nplot.userparams
	while or(nplot.indices(ii)<1,nplot.indices(ii)>info.numparam)
		tmp = input(['which parameter in the ' num2str(ii) 'th position?']);
		while isempty(tmp)
			tmp = input('');
		end
		nplot.indices(ii) = tmp;
	end
end

end