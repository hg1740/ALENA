function Optim = store_envelope_internal_loads(Optim,Parts)

for i = 1:length(Parts)
    
    [Optim.(Parts{i}).Fx_Env(:,1),Optim.(Parts{i}).Fx_Case(:,1)] = max(abs(Optim.(Parts{i}).Fx),[],2); % Axial
    [Optim.(Parts{i}).Fy_Env(:,1),Optim.(Parts{i}).Fy_Case(:,1)] = max(abs(Optim.(Parts{i}).Fy),[],2); % Horizontal
    [Optim.(Parts{i}).Fz_Env(:,1),Optim.(Parts{i}).Fz_Case(:,1)] = max(abs(Optim.(Parts{i}).Fz),[],2); % Vertical
    [Optim.(Parts{i}).Mx_Env(:,1),Optim.(Parts{i}).Mx_Case(:,1)] = max(abs(Optim.(Parts{i}).Mx),[],2); % Torque
    [Optim.(Parts{i}).My_Env(:,1),Optim.(Parts{i}).My_Case(:,1)] = max(abs(Optim.(Parts{i}).My),[],2); % Out plane Bending Moment
    [Optim.(Parts{i}).Mz_Env(:,1),Optim.(Parts{i}).Mz_Case(:,1)] = max(abs(Optim.(Parts{i}).Mz),[],2); % In plane Bending Moment
    
end

end