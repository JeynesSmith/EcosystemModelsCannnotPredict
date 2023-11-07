function Make_TIFF(FileName,Dimensions,Resolution)

% This function creates a TIFF from the current figure
% 
% function Make_TIFF(FileName,Dimensions)
% 
% FileName = the TIFF name and full path
% Dimensions = [X_L Y_L Width Height]
% 


if nargin == 2
   Resolution = '-r200';
end

set(gcf, 'paperunits', 'centimeters', 'paperposition', Dimensions)
set(gcf, 'renderer', 'painters')
print(FileName,'-dtiff',Resolution)
