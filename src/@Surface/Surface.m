classdef Surface
    %SURFACE defines a surface species that can be equilibrated with a
    %solution object. Each surface object should contain only one surface
    %(with or without different types)
    
    properties
        name(1,1) string
        number(1,1) double {mustBeNonnegative, mustBeInteger}
        mass(1,1) double
        scm(1,1) string
        site_density(1,:) double
        surface_master_species(:, 1) string
        surface_species_reactions(:, 1) string
        log_k(:,1) double
        dh(:,1) double
        specific_surface_area(1,1) double {mustBeNonnegative}
        sites_units(1,1) string   % absolute or density
        edl_model(1,1) string
        edl_thickness {mustBeScalarOrEmpty}
        only_counter_ions(1,1) logical
        capacitances(1,2) double
        cd_music_coeffs(:,:) double
    end
    
    methods
        function obj = Surface()
            %SURFACE Construct an instance of this class
            obj.name = "surface";
            obj.number = 1;
        end
        
        function [surface_string, surface_master_string, surface_species_string] = phreeqc_string(obj)
            % phreeqc_string returns a string in phreeqc format based on
            % the defined surface object. It does not run without mixing
            % with a solution object
            % There might be bugs in the conversion and more tests need to
            % be written
            % TBD: currently, surfaces that are linked to and equilibrium
            % or kinetic phases must be defined manually; they are only
            % important is a substantial amount of the surface dissolves or
            % precipitate. For, e.g. flow in a chalk reservoir, it is not
            % significant
            surface_master_string = "SURFACE_MASTER_SPECIES \n";
            for i = 1:length(obj.surface_master_species)
                surface_master_string = strjoin([surface_master_string obj.surface_master_species(i) "\n"]);
            end

            surface_species_string = "SURFACE_SPECIES \n";
            for i = 1:length(obj.surface_species_reactions)
                surface_species_string = strjoin([surface_species_string obj.surface_species_reactions(i) "\n log_k" num2str(obj.log_k(i)) "\n delta_h" num2str(obj.dh(i)) "\n"]);
                if strcmpi(obj.scm, 'cd_music')
                    surface_species_string = strjoin([surface_species_string "-cd_music" num2str(obj.cd_music_coeffs(i,:)) "\n"]);
                end
            end
            
            surface_string = ["SURFACE" num2str(obj.number) obj.name "\n"];
            if strcmpi(obj.scm, 'cd_music')
                surface_string = strjoin([surface_string "-cd_music \n"]);
            end
            if strcmpi(obj.sites_units, 'absolute')
                surface_string = strjoin([surface_string "-sites_units absolute \n"]);
            else
                surface_string = strjoin([surface_string "-sites_units density \n"]);
            end
            ms = strsplit(obj.surface_master_species(1), ' ');
            surface_string = strjoin([surface_string ms(1) num2str(obj.site_density(1)) num2str(obj.specific_surface_area) num2str(obj.mass) "\n"]);
            for i = 2:length(obj.surface_master_species)
                ms = strsplit(obj.surface_master_species(i), ' ');
                surface_string = strjoin([surface_string ms(1) num2str(obj.site_density(i)) "\n"]);
            end
            
            if strcmpi(obj.scm, 'cd_music')
                surface_string = strjoin([surface_string "-capacitances " num2str(obj.capacitances) "\n"]);
            end

            if strcmpi(obj.edl_model, 'diffuse_layer')
                surface_string = strjoin([surface_string "-diffuse_layer " num2str(obj.edl_thickness) "\n"]);
            elseif strcmpi(obj.edl_model, 'Donnan')
                surface_string = strjoin([surface_string "-Donnan " num2str(obj.edl_thickness) "\n"]);
            elseif strcmpi(obj.edl_model, 'no_edl')
                surface_string = strjoin([surface_string "-no_edl \n"]);
            end

            if obj.only_counter_ions
                surface_string = strjoin([surface_string "-only_counter_ions true \n"]);
            end

        end
    end
    
    methods (Static)
        function obj = calcite_surface()
            obj = Surface();
            % TBD
        end
        
        function obj = clay_surface()
            obj = Surface();
            % TBD
        end

        function [surface_string, surface_master_string, surface_species_string] = combine_surfaces(surf1, surf2, new_name, new_number)
            % mixes two surface definitions and creates a phreeqc string
            % that contains both surfaces as a single SURFACE block
            % NOTE: only surfaces with the same model can be combined
            if strcmpi(surf1.scm, surf2.scm)
                % TBD
                [st1, sms1, sss1] = surf1.phreeqc_string();
                [~, sms2, sss2] = surf2.phreeqc_string();
                sms1 = split(sms1, '\n');
                sms2 = split(sms2, '\n');
                surface_master_string = strjoin([sms1, sms2(2:end)], '\n');

                sss1 = split(sss1, '\n');
                sss2 = split(sss2, '\n');
                surface_species_string = strjoin([sss1, sss2(2:end)], '\n');
                
                st1 = split(st1, '\n');
                st1(1) = ["SURFACE" num2str(new_number) new_name "\n"];
                surface_string = strjoin(st1, '\n');
                ms = strsplit(surf2.surface_master_species(1), ' ');
                surface_string = strjoin([surface_string ms(1) num2str(surf2.site_density(1)) num2str(surf2.specific_surface_area) num2str(surf2.mass) "\n"]);
                for i = 2:length(surf2.surface_master_species)
                    ms = strsplit(surf2.surface_master_species(i), ' ');
                    surface_string = strjoin([surface_string ms(1) num2str(surf2.site_density(i)) "\n"]);
                end
            else
                error('Input surfaces use different models and cannot be combined!')
            end
        end

        function obj = read_json(surf_struct)
            % read_json creates a Solution object from a decoded JSON
            % string
            % input:
            %       sol: decoded JSON string to a Matlab structure
            obj = Surface();

            if isfield(surf_struct, 'Name')
                obj.name = surf_struct.Name;
            end

            if isfield(surf_struct, 'Number')
                obj.number = surf_struct.Number;
            end

            if isfield(surf_struct, 'Mass')
                obj.mass = surf_struct.Mass;
            end

            if isfield(surf_struct, 'SCM')
                obj.scm = surf_struct.SCM;
            end

            if isfield(surf_struct, 'SpecificArea')
                obj.specific_surface_area = surf_struct.SpecificArea;
            end

            if isfield(surf_struct, 'MasterSpecies')
                f_names = fieldnames(surf_struct.MasterSpecies); % get the list of components
                s = string(zeros(1, length(f_names)));
                site_dens = zeros(1, length(f_names));
                ms = cellfun(@(x)strjoin(string(getfield(surf_struct.MasterSpecies, x, 'Species'))), fieldnames(surf_struct.MasterSpecies), 'UniformOutput', false);
                sd = cellfun(@(x)getfield(surf_struct.MasterSpecies, x, 'SiteDensity'), fieldnames(surf_struct.MasterSpecies), 'UniformOutput', false);
                s(:) = ms(:);
                site_dens(:) = cell2mat(sd(:));
                obj.surface_master_species = s;
                obj.site_density = site_dens;
                % TBD site_density
            end

            if isfield(surf_struct, 'Reactions')
                f_names = fieldnames(surf_struct.Reactions); % get the list of components
                s = string(zeros(1, length(f_names)));
                dh_val = zeros(1, length(f_names));
                k = zeros(1, length(f_names));
                cdm = zeros(1, length(f_names));
                ms = cellfun(@(x)string(getfield(surf_struct.Reactions, x, 'Reaction')), f_names, 'UniformOutput', false);
                k_cell = cellfun(@(x)getfield(surf_struct.Reactions, x, 'logK'), f_names, 'UniformOutput', false);
                k(:) = cell2mat(k_cell(:));
                try
                    dh_cell = cellfun(@(x)getfield(surf_struct.Reactions, x, 'DeltaH'), f_names, 'UniformOutput', false);
                    dh_val(:) = cell2mat(dh_cell(:));
                catch
                    warning('Values for the enthalpy of surface complexation reactions are not specified; assumed zero.')
                end
                try
                    if strcmpi(obj.scm, 'cd_music')
                        cdm_cell = cellfun(@(x)getfield(surf_struct.Reactions, x, 'cd_music')', f_names, 'UniformOutput', false);
                        cdm = cell2mat(cdm_cell);
                    end
                catch
                    warning('CD MUSIC model coefficients are not assigned.')
                end
                s(:) = ms(:);
                obj.surface_species_reactions = s;
                obj.log_k = k;
                obj.dh = dh_val;
                obj.cd_music_coeffs = cdm;
            end

            if isfield(surf_struct, 'SitesUnits')
                obj.sites_units = surf_struct.SitesUnits;
            end

            if isfield(surf_struct, 'EDL')
                obj.edl_model = surf_struct.EDL;
            end

            if isfield(surf_struct, 'EDL_thickness')
                obj.edl_thickness = surf_struct.EDL_thickness;
            else
                obj.edl_thickness = 1e-8;
            end

            if isfield(surf_struct, 'only_counter_ions')
                obj.only_counter_ions = surf_struct.only_counter_ions;
            else
                obj.only_counter_ions = false;
            end

            if isfield(surf_struct, 'Capacitances')
                obj.capacitances = surf_struct.Capacitances;
            end
        end
    end
end

