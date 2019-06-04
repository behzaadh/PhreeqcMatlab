function phreeqc_rm = PhreeqcSingleCell(input_file, data_base)
% PHREEQCSINGLECELL Creates a phreeqcrm intance that works on a single reaction cell
%   input_file: the full name or address to the phreeqc input file
% in the input file, all the phase, surface, exchange, etc must be numbered
% as number 1 for this function to work properly. See the example input
% file. The input file must be clean at the moment. No commenting out the
% lines, although I do a bit of clean up in the input file.

phreeqc_rm = PhreeqcRM(1, 1); % one cell, one thread
phreeqc_rm = phreeqc_rm.RM_Create(); % create a PhreeqcRM instance
status = phreeqc_rm.RM_SetComponentH2O(false);
status = phreeqc_rm.RM_SetUnitsSolution(2);           % 1, mg/L; 2, mol/L; 3, kg/kgs
status = phreeqc_rm.RM_SetUnitsPPassemblage(1);       % 0, mol/L cell; 1, mol/L water; 2 mol/L rock
status = phreeqc_rm.RM_SetUnitsExchange(1);           % 0, mol/L cell; 1, mol/L water; 2 mol/L rock
status = phreeqc_rm.RM_SetUnitsSurface(1);            % 0, mol/L cell; 1, mol/L water; 2 mol/L rock
status = phreeqc_rm.RM_SetUnitsGasPhase(1);           % 0, mol/L cell; 1, mol/L water; 2 mol/L rock
status = phreeqc_rm.RM_SetUnitsSSassemblage(1);       % 0, mol/L cell; 1, mol/L water; 2 mol/L rock
status = phreeqc_rm.RM_SetUnitsKinetics(1);           % 0, mol/L cell; 1, mol/L water; 2 mol/L rock
phreeqc_rm.RM_UseSolutionDensityVolume(true);

status = phreeqc_rm.RM_SetPorosity(1.0);             % If pororosity changes due to compressibility
status = phreeqc_rm.RM_SetSaturation(1.0);           % If saturation changes
    
status = phreeqc_rm.RM_LoadDatabase(database_file(data_base)); % load the database
status = phreeqc_rm.RM_RunFile(true, true, true, input_file); % run the input file

ncomps = phreeqc_rm.RM_FindComponents();
comp_name = phreeqc_rm.GetComponents();

ic1 = -1*ones(7, 1);
ic2 = -1*ones(7, 1);
f1 = ones(7, 1);

% look for the keywords in the inputfile
C = ReadPhreeqcFile(input_file); % read and clean the input file

if any(contains(C, 'SELECTED_OUTPUT')) % Surface 1
    status = phreeqc_rm.RM_SetSelectedOutputOn(true);
end

if ~any(contains(C, 'SOLUTION'))
    error('PhreeqcMatlab: SOLUTION 1 must be defined in the input file.');
end

ic1(1) = 1;              % Solution 1

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

status = phreeqc_rm.RM_InitialPhreeqc2Module(ic1, ic2, f1);
status = phreeqc_rm.RM_RunCells();

end

