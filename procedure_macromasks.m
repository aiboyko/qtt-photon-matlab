%% hardcoded macromasks (needed for multiscale geometries)
% macromask0 = [   1 1 1 0 1 1 1 1;...
%                  1 1 1 0 1 1 1 1;...
%                  1 1 1 0 1 1 1 1;...
%                  1 1 1 0 1 1 1 1;...
%                  0 0 0 0 0 0 0 0;...
%                  1 1 1 0 1 1 1 1;...
%                  1 1 1 0 1 1 1 1;...
%                  1 1 1 0 1 1 1 1];
%              
%                  
% macromask = [0 0 1 1 0 0 1 1;...
%              0 0 1 1 0 0 1 1;...
%              1 1 0 0 1 1 0 0;...
%              1 1 0 0 1 1 0 0;...
%              0 0 1 1 0 0 1 1;...
%              0 0 1 1 0 0 1 1;...
%              1 1 0 0 1 1 0 0;...
%              1 1 0 0 1 1 0 0];
%          
% macromask1 = [0 1 0 0 0 0 0 0;...
%               0 1 1 1 1 1 1 1;...
%               0 1 1 0 0 0 1 0;...
%               0 1 0 0 1 0 1 0;...
%               0 1 0 1 0 0 1 0;...
%               0 1 0 0 0 1 1 0;...
%               1 1 1 1 1 1 1 0;...
%               0 0 0 0 0 0 1 0];
%           
% macromask1 = [0 1 0 0 0 0 1 0;...
%               1 1 0 0 0 0 1 1;...
%               0 0 0 0 0 1 0 0;...
%               0 0 0 0 1 0 0 0;...
%               0 0 0 1 0 0 0 0;...
%               0 0 1 0 0 0 0 0;...
%               1 1 0 0 0 0 0 0;...
%               0 1 0 0 0 0 0 0];           
%          
% macromask=kron(macromask0,macromask1);          