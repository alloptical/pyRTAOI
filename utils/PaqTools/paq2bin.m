[data, names, units, rate,filename]=paq2lab;
FileToWrite=strrep(filename,'.paq','.bin');
fid=fopen(FileToWrite,'w+');
fwrite(fid,data','float32');
fclose(fid);