cd MATLAB_files/
files=dir('*.mat');
cd ..
for i=1:numel(files)
    model_name=files(i).name;
    completed=mat2tf(model_name);
end