classdef solution
    %SOLUTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        number
        unit
        components
        concentrations
        charge_balance_component
        charge_balance
        density
        pH
        pe
        alkalinity
        alkalinity_component
    end
    
    methods
        function obj = solution(primary_species, concentration)
            %SOLUTION Construct an instance of this class
            % Input arguments:
            % primary_species: an array of strings
            % concentration: an array of concentration values
            %   Detailed explanation goes here
            obj.components = primary_species;
            obj.concentrations = concentration;
        end
        
        function obj = seawater(obj)
            obj.name = "Seawater";
            obj.number = 1;
            obj.unit = "ppm";
            obj.components = ["Ca", "Mg", "Na", "K", "Si", "Cl", "S(6)"];
            obj.concentrations = [412.3, 1291.8, 10768.0, 399.1, 4.28, 19353.0, 2712.0];
            obj.charge_balance_component = [];
            obj.charge_balance = false;
            obj.density = 1.0253; % kg/l
            obj.pH = 8.22;
            obj.pe = 8.451;
            obj.alkalinity = 141.682;
            obj.alkalinity_component = "HCO3";
%                 units   ppm
%                 pH      8.22
%                 pe      8.451
%                 density 1.023
%                 temp    25.0
%                 Ca              412.3
%                 Mg              1291.8
%                 Na              10768.0
%                 K               399.1
%                 Si              4.28
%                 Cl              19353.0
%                 Alkalinity      141.682 as HCO3
%                 S(6)            2712.0
        end
        
        function solution_string = phreeqc_solution(obj)
            %phreeqc_solution returns a string of phreeqc format for the
            %string object
            outputArg = obj.Property1 + inputArg;
        end
    end
end

