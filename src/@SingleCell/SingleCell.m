classdef SingleCell
    %SINGLECELL is a class that contains a mixture of different components.
    % it contains at least an aqueous solution. It can contain solid and
    % gas phases too, which are defined as PhreeqcMatlab objects. It also
    % has fields for temperature and pressure.
    % The reaction rate field can also be specified to run the cell for a
    % certain specified time. Running a cell does not update its fields
    % (i.e. equilibrate aqueous solutions with other specified phases, etc)
    % but one can get a copy of an aqueous cell after running it with
    % PhreeqcRM. The copy is saved as a PhreeqcCell class.
    
    properties
        temperature
        pressure
        aqueous_solution
        equilibrium_phase
        surface_specie
        exchanger
        reaction_rate
        gas_phase
    end
    
    methods
        function obj = SingleCell(aq_solution, varargin)
            %SINGLECELL constructs a single cell that contains an aqueous
            %liquid and other phases, e.g. soluble minerals, surfaces,
            %exchangers, gas phases, etc.
            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

