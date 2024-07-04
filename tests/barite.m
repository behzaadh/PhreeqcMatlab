
C_components = [55.5253251007785,0,0,0,0.999999999999998,0.0100000000000002,0.0999999999999998,0.0999999999999998,0.0999999999999998,0.00100000000000019;...
    55.5253251007785,0,0,0,0.999999999999962,0.0100000000000034,0.0999999999999966,0.0999999999999966,0.0999999999999966,0.00100000000000373;...
    55.5253251007785,0,0,0,0.999999999999188,0.0100000000000731,0.0999999999999269,0.0999999999999269,0.0999999999999269,0.00100000000008041;...
    55.5253251007785,0,0,0,0.999999999981096,0.0100000000017014,0.0999999999982986,0.0999999999982986,0.0999999999982986,0.00100000000187150;...
    55.5253251007785,0,0,0,0.999999999522609,0.0100000000429656,0.0999999999570344,0.0999999999570344,0.0999999999570344,0.00100000004726217;...
    55.5253251007785,0,0,0,0.999999986852127,0.0100000011833204,0.0999999988166796,0.0999999988166796,0.0999999988166796,0.00100000130165246;...
    55.5253251007785,0,0,0,0.999999602677431,0.0100000357593888,0.0999999642406112,0.0999999642406112,0.0999999642406112,0.00100003933532771;...
    55.5253251007785,0,0,0,0.999986727917285,0.0100011944993894,0.0999988055006106,0.0999988055006106,0.0999988055006106,0.00100131394932830;...
    55.5253251007785,0,0,0,0.999505492364312,0.0100445061322733,0.0999554938677268,0.0999554938677268,0.0999554938677268,0.00104895674550058;...
    55.5253251007785,0,0,0,0.979214498573520,0.0118707138355216,0.0981292861644784,0.0981292861644784,0.0981292861644784,0.00305778521907373];

T = [76.6666666666667;76.6666666666667;76.6666666666667;76.6666666666667;76.6666666666667;76.6666666666667;76.6666666666667;76.6665988498264;76.6608513726128;76.1097971598307];
threads = 1;
database_name = '../database/phreeqc.dat';
pgcfile = 'barite.pqi';
nCells = size(C_components, 1);
phreeqc_rm = PhreeqcRM(nCells, threads);
phreeqc_rm = phreeqc_rm.RM_Create();
phreeqc_rm.RM_UseSolutionDensityVolume(true);

status = phreeqc_rm.RM_LoadDatabase(database_name);

status = phreeqc_rm.RM_RunFile(true, true, true, pgcfile);

status = phreeqc_rm.RM_SetSelectedOutputOn(true);

ic1 = -1*ones(nCells*7, 1);
ic2 = -1*ones(nCells*7, 1);
f1 = ones(nCells*7, 1);
for i = 1:nCells
    ic1(i) = 1;              % Solution 1
    ic1(1*nCells + i) = +1;    % Equilibrium phases
    ic1(2*nCells + i) = -1;    % Exchange 1
    ic1(3*nCells + i) = -1;    % Surface none
    ic1(4*nCells + i) = -1;    % Gas phase none
    ic1(5*nCells + i) = -1;    % Solid solutions none
    ic1(6*nCells + i) = -1;    % Kinetics none
end
status = phreeqc_rm.RM_InitialPhreeqc2Module(ic1, ic2, f1);

nPPAssemblage = phreeqc_rm.RM_GetPPAssemblageCount();
PPAssemblageComps = phreeqc_rm.GetPPAssemblageComps();
moles= [1;2;3;4;1;6;7;8;9;10;1;2;3;4;1;6;7;8;9;10;1;2;3;4;1;6;7;8;9;10;1;2;3;4;1;6;7;8;9;10;1;2;3;4;1;6;7;8;9;10];
status = phreeqc_rm.RM_SetPPAssemblageMoles(moles);
si = moles/10;
status = phreeqc_rm.RM_SetPPAssemblageSI(si);


nComps = phreeqc_rm.RM_FindComponents();
names_comps = phreeqc_rm.GetComponents()';

nGasComps = phreeqc_rm.RM_GetGasComponentsCount();
names_gas_comps = phreeqc_rm.GetGasComponentsNames()';

%% section added by Ali starting 24-1-2024
status = phreeqc_rm.RM_SetSpeciesSaveOn(1);

PORO = 0.3 * ones(1,nCells); 
status = phreeqc_rm.RM_SetPorosity(PORO);

Sw = 1 * ones(1,nCells);
status = phreeqc_rm.RM_SetSaturation(Sw);

v = 1 * ones(1,nCells); % liters
status = phreeqc_rm.RM_SetRepresentativeVolume(v);

status = phreeqc_rm.RM_SetUnitsSolution(2); %Options are 1, mg/L; 2 mol/L; or 3, mass fraction, kg/kgs.

C_components = reshape(C_components, 1, nCells*nComps);
    
status = phreeqc_rm.RM_SetConcentrations(C_components);
status = phreeqc_rm.RM_SetTemperature(T);
%% Run
status = phreeqc_rm.RM_RunCells();

heading_out = phreeqc_rm.GetSelectedOutputHeadings(1);
selected_out = phreeqc_rm.GetSelectedOutput(1);

c_out_comps = phreeqc_rm.GetConcentrations();

names_species = phreeqc_rm.GetSpeciesNames();
c_out_species = phreeqc_rm.GetSpeciesConcentrations();
    
status_destroy = phreeqc_rm.RM_Destroy();

C_components = reshape(c_out_comps, nCells, nComps);