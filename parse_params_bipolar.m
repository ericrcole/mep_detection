function [inds] = parse_params_bipolar(mode, param_struct, ep_data, amp_type, switch_lead, lead_type)
    %loop through stim param strings and compile into matrix for plotting
    %based on desired mode (monopolar segmented, bipolar, etc)
    %
    %mode = which montage to use
    %param_struct = data struct with stim parameter strings, 
    
    switch mode
        case 'monopolar_segmented'
            %stim_dict = containers.Map({-1,1,2,3,4,21,22,23,31,32,33},{-1 1 nan nan 8 2 3 4 5 6 7});
            stim_dict = containers.Map({'C','1','2A','2B','2C','3A','3B','3C','4','2','3'},{-1,1,2,3,4,5,6,7,8,nan,nan});
            %% Checking amp type
            if strcmp(amp_type,'any')
                inds = nan(8,1); % creates an 8 by 1 array of NaN
            elseif isvector(amp_type)
                inds = nan(8,length(amp_type)); % creates an 8 by however many amps given array of NaN (usually 8 * 3)
                amps = amp_type;
            else % fallthrough case that will calculate directly from the structure.
                amps = abs(unique(param_struct.amplitudes));
                inds = nan(8,length(amps));
            end
            
            %% Looks through parameter strings
            for k = 1:length(param_struct.param_strings)

                amp = abs(param_struct.amplitudes(k));
                if (iscell(param_struct.cathodes{k}) && length(param_struct.cathodes{k})>1) || (iscell(param_struct.anodes{k}) && length(param_struct.anodes{k})>1)
                    continue
                end
                if iscell(param_struct.cathodes{k})
                    cathode = stim_dict(param_struct.cathodes{k}{1});
                else
                    cathode = stim_dict(param_struct.cathodes{k});
                end
                if iscell(param_struct.anodes{k})
                    anode = stim_dict(param_struct.anodes{k}{1});
                else
                    anode = stim_dict(param_struct.anodes{k});
                end
                % cathode and anode should be mapped to a value according
                % to stim_dict. This value can serve as a direct index for
                % the final matrix.
                
                if (strcmp(amp_type,'any') && anode == -1 && ~isnan(cathode)) %% If chest-implanted monitor was included, which won't be in bipolar case.
                    inds(cathode) = ep_data(k);
                elseif (anode == -1) && ~isnan(cathode)
                    ampind = find(amps == amp,1);
                    inds(cathode,ampind) = ep_data(k);
                end

            end
            if switch_lead && contains(lead_type, 'ABB')
                if size(inds,2) > 1
                    inds = inds([1 2 4 3 5 7 6 8],:);
                else
                    inds = inds([1 2 4 3 5 7 6 8]);
                end
            end
            
        case 'monopolar'
            stim_dict = containers.Map({'C','1','2A','2B','2C','3A','3B','3C','4','2','3'},{-1,1,nan,nan,nan,nan,nan,nan,4,2,3});
            %stim_dict = containers.Map({-1,1,2,3,4,21,22,23,31,32,33},{nan 1 2 3 4 nan nan nan nan nan nan});
            if strcmp(amp_type,'any')
                inds = zeros(4,1);
            elseif isvector(amp_type)
                amps = amp_type;
                inds = nan(4,length(amps));
            else
                amps = unique(param_struct.amplitudes);
                inds = zeros(4,length(amps));
            end
            
            for k = 1:length(param_struct.param_strings)
                    
                amp = param_struct.amplitudes(k);
%                 if (iscell(param_struct.cathodes{k}) && length(param_struct.cathodes{k})>1) || (iscell(param_struct.anodes{k}) && length(param_struct.anodes{k})>1)
%                     continue
%                 end
                if iscell(param_struct.cathodes{k}) && length(param_struct.cathodes{k})==3
                    cathode = stim_dict(param_struct.cathodes{k}{1}(1));
                elseif iscell(param_struct.cathodes{k}) && length(param_struct.cathodes{k})==1
                    cathode = stim_dict(param_struct.cathodes{k}{1});
                elseif iscell(param_struct.cathodes{k}) && length(param_struct.cathodes{k})>1
                    cathode = nan;
                else
                    cathode = stim_dict(param_struct.cathodes{k});
                end
                
                if iscell(param_struct.anodes{k}) && length(param_struct.anodes{k})==3
                    anode = stim_dict(param_struct.cathodes{k}{1}(1));
                elseif iscell(param_struct.anodes{k}) && length(param_struct.anodes{k})==1
                    anode = stim_dict(param_struct.anodes{k}{1});
                elseif iscell(param_struct.anodes{k}) && length(param_struct.anodes{k})>1
                    anode = nan;
                else
                    anode = stim_dict(param_struct.anodes{k});
                end
                
                if (strcmp(amp_type,'any') && anode == -1 && ~isnan(cathode))
                    inds(cathode) = ep_data(k);
                
                elseif (anode == -1) && ~isnan(cathode)
                    ampind = find(amps == amp,1);
                    inds(cathode,ampind) = ep_data(k);
                end
                
            end
        case 'bipolar_segmented'
            %% Create mapping object
            if strcmp(lead_type, 'BSC')
                stim_dict = containers.Map({'C','1','2','3','4','5','6','7','8','L2','L3'},{-1,1,2,3,4,5,6,7,8,nan,nan});
            else
                stim_dict = containers.Map({'C','1','2A','2B','2C','3A','3B','3C','4','2','3'},{-1,1,2,3,4,5,6,7,8,nan,nan});
            end
            %% Check AMP_TYPE
            if strcmp(amp_type,'any')
                inds = nan(8,8,1); % creates an 8 by 1 array of NaN
            elseif isvector(amp_type)
                inds = nan(8,8,length(amp_type)); % creates an 8 by however many amps given array of NaN (usually 8 * 3)
                amps = amp_type;
            else % fallthrough case that will calculate directly from the structure.
                amps = abs(unique(param_struct.amplitudes));
                inds = nan(8,8,length(amps));
            end
            for k = 1:length(param_struct.param_strings)

                amp = abs(param_struct.amplitudes(k));
                if (iscell(param_struct.cathodes{k}) && length(param_struct.cathodes{k})>1) || (iscell(param_struct.anodes{k}) && length(param_struct.anodes{k})>1)
                    continue
                end
                if iscell(param_struct.cathodes{k})
                    cathode = stim_dict(param_struct.cathodes{k}{1});
                else
                    cathode = stim_dict(param_struct.cathodes{k});
                end
                if iscell(param_struct.anodes{k})
                    anode = stim_dict(param_struct.anodes{k}{1});
                else
                    anode = stim_dict(param_struct.anodes{k});
                end
                % cathode and anode should be mapped to a value according
                % to stim_dict. This value can serve as a direct index for
                % the final matrix.
                
                if (strcmp(amp_type,'any') && ~isnan(anode) && ~isnan(cathode)) %% If chest-implanted monitor was included, which won't be in bipolar case.
                    inds(cathode, anode, 1) = ep_data(k);
                elseif ~isnan(anode) && ~isnan(cathode) && anode > 0
                    ampind = find(amps == amp,1); % which amp is it?
                    inds(cathode, anode, ampind) = ep_data(k);
                end

            end
        case 'bipolar'
            stim_dict = containers.Map({'C','1','2A','2B','2C','3A','3B','3C','4','2','3'},{nan,1,nan,nan,nan,nan,nan,nan,4,2,3});
            % Improvements for segmented code could include whether or not
            % 2A2B (but not C) was stimulated or 2A-2B-2C was listed
            % instead of 2. Clinical consistency would be ideal.
            %% CHECK AMP_TYPE
            if strcmp(amp_type,'any')
                inds = nan(4,4,1); % creates an 4 by 4 by 1 array of NaN --> cathode * anode * amp count
            elseif isvector(amp_type)
                inds = nan(4,4,length(amp_type)); 
                amps = amp_type;
            else % fallthrough case that will calculate directly from the structure.
                amps = abs(unique(param_struct.amplitudes));
                inds = nan(4,4,length(amps));
            end
            for k = 1:length(param_struct.param_strings)
                    
                amp = param_struct.amplitudes(k);
%                 if (iscell(param_struct.cathodes{k}) && length(param_struct.cathodes{k})>1) || (iscell(param_struct.anodes{k}) && length(param_struct.anodes{k})>1)
%                     continue
%                 end
                if iscell(param_struct.cathodes{k}) && length(param_struct.cathodes{k})==3
                    cathode = stim_dict(param_struct.cathodes{k}{1}(1));
                elseif iscell(param_struct.cathodes{k}) && length(param_struct.cathodes{k})==1
                    cathode = stim_dict(param_struct.cathodes{k}{1});
                elseif iscell(param_struct.cathodes{k}) && length(param_struct.cathodes{k})>1
                    cathode = nan;
                else
                    cathode = stim_dict(param_struct.cathodes{k});
                end
                
                if iscell(param_struct.anodes{k}) && length(param_struct.anodes{k})==3
                    anode = stim_dict(param_struct.cathodes{k}{1}(1));
                elseif iscell(param_struct.anodes{k}) && length(param_struct.anodes{k})==1
                    anode = stim_dict(param_struct.anodes{k}{1});
                elseif iscell(param_struct.anodes{k}) && length(param_struct.anodes{k})>1
                    anode = nan;
                else
                    anode = stim_dict(param_struct.anodes{k});
                end
                

                if (strcmp(amp_type,'any') && ~isnan(anode) && ~isnan(cathode)) %% If chest-implanted monitor was included, which won't be in bipolar case.
                    inds(cathode, anode, 1) = ep_data(k);
                elseif ~isnan(anode) && ~isnan(cathode) && anode > 0
                    ampind = find(amps == amp,1); % which amp is it?
                    inds(cathode, anode, ampind) = ep_data(k);
                    
                elseif ~isnan(anode) && ~isnan(cathode) && anode == -1
                    ampind = find(amps == amp,1); % which amp is it?
                    inds(cathode, cathode, ampind) = ep_data(k);
                end
            end
    end

    
    
end