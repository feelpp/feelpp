% Constrained circle packing example
ab=[500 500]; % rectangle dimensions
R_min=5;      % minimum circle radius
R_max=25;     % maximum circle radius
cnst=true;   
[C,R]=random_circle_packing_rectangle(ab,R_min,R_max,cnst);