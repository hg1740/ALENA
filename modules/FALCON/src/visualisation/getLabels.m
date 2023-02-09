function textout = getLabels(dirn,forceorvel_flag,unitflag,prefix,suffix)

if nargin < 2 || isempty(forceorvel_flag)
    labels = {'Axial Shear','In-plane Shear','Vertical Shear','Torque','Bending','In-plane Bending'};
    units  = {'N','N','N','Nm','Nm','Nm'};
elseif strcmp(forceorvel_flag(1),'f')
    if length(forceorvel_flag) < 2 || strcmp(forceorvel_flag(2),'z')
        labels = {'Axial Shear','In-plane Shear','Vertical Shear','Torque','Bending','In-plane Bending'};
        units  = {'N','N','N','Nm','Nm','Nm'};
    else
        labels = {'Axial Shear','Vertical Shear','In-plane Shear','Torque','In-plane Bending','Bending'};
        units  = {'N','N','N','Nm','Nm','Nm'};
    end
elseif strcmp(forceorvel_flag(1),'v')
    if length(forceorvel_flag) < 2 || strcmp(forceorvel_flag(2),'y')
        labels = {'Lat velocity','Forward velocity','Vert velocity','Pitch','Roll','Yaw'};
        units  = {'ms^{-1}','ms^{-1}','ms^{-1}','rads^{-1}','rads^{-1}','rads^{-1}'};
    else
        labels = {'Forward velocity','Lat velocity','Vert velocity','Roll','Pitch','Yaw'};
        units  = {'ms^{-1}','ms^{-1}','ms^{-1}','rads^{-1}','rads^{-1}','rads^{-1}'};
    end
else
    error('Don''t know how to interpret that!')
end

textout = labels{dirn};

if nargin > 3 && ~isempty(prefix)
    textout = [prefix,' ',textout];
end
if nargin > 4 && ~isempty(suffix)
    textout = [textout,' ',suffix];
end
if nargin > 2 && unitflag == 0
else
    textout = [textout,' (',units{dirn},')'];
end