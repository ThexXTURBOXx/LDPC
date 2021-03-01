% Read all pictures into MATLAB
folder  = '.';
list    = dir(fullfile(folder, 'pic*.m'));
nFile   = length(list);
success = false(1, nFile);
for k = 1:nFile
  % Execute file
  file = list(k).name;
  try
    run(fullfile(folder, file));
    success(k) = true;
  catch
    fprintf('failed: %s\n', file);
  end

  TEMP = reshape(TEMP(1:10000), 100, 100)';
  imwrite(TEMP, strcat(file, '.png'));

  % Clear picture
  clear TEMP;
end

