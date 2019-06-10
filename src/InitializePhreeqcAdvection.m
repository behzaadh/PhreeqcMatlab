function phreeqc_rm = InitializePhreeqcAdvection(phreeqc_rm, input_file)
%INITIALIZEPHREEQCADVECTION Initializes a phreeqc instance with the phreeqc
%input file for inital and boundary conditions

% map transport grid to the reaction cells (1 to 1 since it is only 1D)
nxyz = phreeqc_rm.ncells;
grid2chem = 0:nxyz-1; % o base indexing for C
status = phreeqc_rm.RM_CreateMapping(grid2chem);

workers = true;             % Worker instances do the reaction calculations for transport
initial_phreeqc = true;     % InitialPhreeqc instance accumulates initial and boundary conditions
utility = true;             % Utility instance is available for processing
status = phreeqc_rm.RM_RunFile(workers, initial_phreeqc, utility, 'advect.pqi');
% Clear contents of workers and utility
initial_phreeqc = false;
string_input = 'DELETE; -all';
status = phreeqc_rm.RM_RunString(workers, initial_phreeqc, utility, string_input);
% Determine number of components to transport
ncomps = phreeqc_rm.RM_FindComponents();

% Initial condition
% it is currently not possible to get the total number of solution,
% surface, etc blocks in the phreeqc input file, automatically. Therefore,
% I'm assuming that if a keyword exist in the input file, nxyz block of
% that keyword is defined in the input file. Later, I will think of a
% better method (probably a matlab expression in the transport input file)
ic1 = -1*ones(nxyz, 7);
ic2 = -1*ones(nxyz, 7);
f1 = ones(nxyz, 7);
ic1(:,1) = 1;              % Solution 1
ic1(:,2) = -1;      % Equilibrium phases none
ic1(:,3) = 1;     % Exchange 1
ic1(:,4) = -1;    % Surface none
ic1(:,5) = -1;    % Gas phase none
ic1(:,6) = -1;    % Solid solutions none
ic1(:,7) = -1;    % Kinetics none

status = phreeqc_rm.RM_InitialPhreeqc2Module(ic1, ic2, f1);







C = ReadPhreeqcFile(input_file); % read and clean the input file

if any(contains(C, 'SELECTED_OUTPUT')) % Surface 1
    status = phreeqc_rm.RM_SetSelectedOutputOn(true);
end

if ~any(contains(C, 'SOLUTION'))
    error('PhreeqcMatlab: SOLUTION must be defined in the input file.');
end

ic1(1) = 1:nxyz;              % Solution 1-n

% in phreeqc: EQUILIBRIUM_PHASES is the keyword for the data block. Optionally, EQUILIBRIUM , EQUILIBRIA , PURE_PHASES , PURE .
if any(contains(C, 'EQUILIBRIUM_PHASES')) ||  any(contains(C, 'EQUILIBRIUM')) || any(contains(C, ' EQUILIBRIA')) || any(contains(C, ' PURE_PHASES')) || any(contains(C, ' PURE'))
    ic1(2) = 1;      % Equilibrium phases
end

% Exchange species
if any(contains(C, 'EXCHANGE'))
    ic1(3) = 1;     % Exchange 1
end

% Surface species

if any(contains(C, 'SURFACE')) 
    ic1(4) = 1; % Surface 1
end

% Gas phase
if any(contains(C, 'GAS_PHASE')) % Surface 1
    ic1(5) = 1;    % Gas phase 1
end

% Solid solution
if any(contains(C, 'SOLID_SOLUTION')) % Surface 1
    ic1(6) = 1;    % Solid solutions 1
end

% Kinetics
if any(contains(C, 'KINETICS')) % Surface 1
    ic1(7) = 1;    % Kinetics 1
end

end

