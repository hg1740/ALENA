function beam_model =  GenerateBCs(beam_model,Aircraft_param)
%   - C.Szczyglowski(08/10/2019) This is already performed by the
%   'convertToFE' method of 'FEable'.

Aircraft_param_parts = fieldnames(Aircraft_param);

RBE2Index = 1000;

for i = 1:length(Aircraft_param_parts)
    if isfield(Aircraft_param.(Aircraft_param_parts{i}),'LiftingSurface')
        if Aircraft_param.(Aircraft_param_parts{i}).LiftingSurface== 1
            % If its a lifting surface connect it to the parent
            % If the parent is a lifting surface then generate cut one of the
            % beams and connect it, otherwise connect it rigidly to body
            
            Parent = Aircraft_param.(Aircraft_param_parts{i}).Parent;
            
            if ~isempty(Parent)
                
                if Aircraft_param.(Parent).LiftingSurface == 1
                    error('Currently unable to connect a lifting surface to another lifting surface');
                else
                    
                    if isfield(Aircraft_param.(Aircraft_param_parts{i}).Planform,'CT_y_eta')
                        
                        % Create Rigid connection between the centreline of the part
                        % and the parent
                        [beam_model.RBE2,beam_model.Info,beam_model.PartIDs.(Aircraft_param_parts{i}).RBE2IDs]  = ...
                            RBE2Creator(beam_model.RBE2,beam_model.Info,RBE2Index + i*100,Aircraft_param.(Parent).Planform.NodeIdx + 1,...
                            Aircraft_param.(Aircraft_param_parts{i}).Planform.NodeIdx+1,'246',[]);
                        
                        % Create Rigid connection between the root of the part
                        % and the parent
                        
                        [beam_model.RBE2,beam_model.Info,beam_model.PartIDs.(Aircraft_param_parts{i}).RBE2IDs]  = ...
                            RBE2Creator(beam_model.RBE2,beam_model.Info,RBE2Index + i*100 + 1,Aircraft_param.(Parent).Planform.NodeIdx + 1,...
                            Aircraft_param.(Aircraft_param_parts{i}).Planform.CTNode,'135',[]);
                        
                    else
                        
                        [beam_model.RBE2,beam_model.Info,beam_model.PartIDs.(Aircraft_param_parts{i}).RBE2IDs]  = ...
                            RBE2Creator(beam_model.RBE2,beam_model.Info,RBE2Index + i*100,Aircraft_param.(Parent).Planform.NodeIdx + 1,...
                            Aircraft_param.(Aircraft_param_parts{i}).Planform.NodeIdx+1,'123456',[]);
                    end
                end
            end
            
        else
            
            Parent = Aircraft_param.(Aircraft_param_parts{i}).Parent;
            
            % Place the bluff body node as half way through the body
            Mid_node = Aircraft_param.(Aircraft_param_parts{i}).Planform.ref_point + [Aircraft_param.(Aircraft_param_parts{i}).Planform.Length*Aircraft_param.(Aircraft_param_parts{i}).Planform.y_eta(end)/2,0,0];
            
            if isempty(Parent)
                
                [beam_model.Node,beam_model.Info,beam_model.PartIDs.(Aircraft_param_parts{i}).BeamNodes]  = ...
                    NodeCreator(beam_model.Node,Mid_node(1),Mid_node(2),Mid_node(3),beam_model.Info, Aircraft_param.(Aircraft_param_parts{i}).Planform.NodeIdx);
            end
            
        end
        
        % If the part has no parents, then it is the parent of all!
        if isempty(Parent)
            
            [beam_model.SPC,beam_model.Info]   = ...
                SPCCreator(beam_model.SPC,beam_model.Info, Aircraft_param.(Aircraft_param_parts{i}).Planform.NodeIdx + 1,...
                123456, Aircraft_param.(Aircraft_param_parts{i}).Planform.NodeIdx + 1);
            
            beam_model.Param.SPC = Aircraft_param.(Aircraft_param_parts{i}).Planform.NodeIdx + 1;
        end
    end
end

end