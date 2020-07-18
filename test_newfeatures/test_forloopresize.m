psy=1:100;
i=1;
while i<numel(psy)
    i=i+1
    psy=psy(1:numel(psy)-1);
end
%so it DOESNT work that way (for is pre-defined in the beginning)
% while loop, however, does re-check