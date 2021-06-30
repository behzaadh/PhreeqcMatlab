classdef surface
    %SURFACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        number
        binding_site
        surface_master_species
        surface_species_reactions
        log_k
        dh
        specific_surface_area
        eq_solution
        sites_units   % absolute or density
        edl_layer
    end
    
    methods
        function obj = surface(binding_site,surface_master_species)
            %SURFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.binding_site = binding_site;
            obj.surface_master_species = surface_master_species;
        end
        
        function surface_string = phreeqc_string(obj)
           surface_string = ["SURFACE" num2str(obj.number) obj.name "\n"];
        end
    end
    
    methods (Static)
        function obj = calcite_surface()
            obj = surface();
        end
        
        function obj = clay_surface()
            obj = surface();
        end
    end
end

