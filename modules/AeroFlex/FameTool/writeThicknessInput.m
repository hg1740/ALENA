function writeThicknessInput(homeFolder,obj)

% Initialise the file where the data will be stored
fp = fopen(fullfile(homeFolder,'\ThicknessInput.m'), 'w');

fprintf(fp,'%% Initialise the thickness\n');

% Write the skin thickness distribution

[~,idx] = unique(obj.Fame.Results.Thickness.eta);
eta_ = regexprep(mat2str(obj.Fame.Results.Thickness.eta(idx)),'\s*',',');
thi_ = regexprep(mat2str(obj.Fame.Results.Thickness.Ush(idx)/1000),'\s*',',');
max_ = sprintf('%s',num2str(0.15));
min_ = sprintf('%s',num2str(0.003));

fprintf(fp,'param.Wing.Box.skin.eta       = %s;\n',eta_);
fprintf(fp,'param.Wing.Box.skin.thickness = %s;\n',thi_);
fprintf(fp,'param.Wing.Box.skin.max       = %s;\n',max_);
fprintf(fp,'param.Wing.Box.skin.min       = %s;\n',min_);
fprintf(fp,'\n');

% Write the spar thickness distribution
[~,idx] = unique(obj.Fame.Results.Thickness.eta);
eta_ = regexprep(mat2str(obj.Fame.Results.Thickness.eta(idx)),'\s*',',');
thi_ = regexprep(mat2str(obj.Fame.Results.Thickness.Fsp(idx)/1000),'\s*',',');
max_ = sprintf('%s',num2str(0.15));
min_ = sprintf('%s',num2str(0.003));
fprintf(fp,'param.Wing.Box.spar.eta       = %s;\n',eta_);
fprintf(fp,'param.Wing.Box.spar.thickness = %s;\n',thi_);
fprintf(fp,'param.Wing.Box.spar.max       = %s;\n',max_);
fprintf(fp,'param.Wing.Box.spar.min       = %s;\n',min_);
fprintf(fp,'\n');

% Write the stringer area distribution
eta_ = regexprep(mat2str([0;1]),'\s*',',');
thi_ = regexprep(mat2str([0.0001;0.0001]),'\s*',',');
max_ = sprintf('%s',num2str(0.05));
min_ = sprintf('%s',num2str(0.00005));
fprintf(fp,'param.Wing.Box.stringer.eta   = %s;\n',eta_);
fprintf(fp,'param.Wing.Box.stringer.area  = %s;\n',thi_);
fprintf(fp,'param.Wing.Box.stringer.max   = %s;\n',max_);
fprintf(fp,'param.Wing.Box.stringer.min   = %s;\n',min_);
fprintf(fp,'\n');
fclose(fp);

end
