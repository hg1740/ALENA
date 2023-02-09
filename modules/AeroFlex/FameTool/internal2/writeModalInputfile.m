function writeModalInputfile(homeFolder)

fp = fopen(fullfile(homeFolder,'ModalInputFile.inc'), 'w');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'$------- MODE PROPERTIES -------------------------------------------------\n');
fprintf(fp,'$-------2-------3-------4-------5-------6-------7-------8-------9-------10\n');
fprintf(fp,'EIGR    1               0       999999          70      \n');
fprintf(fp,'        MASS \n');
fclose(fp);
end
