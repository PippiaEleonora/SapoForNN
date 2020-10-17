function [completed] = mat2tf(model_name)
%mat2tf Translates a MATLAB structure into TF format
%   The NN model might a set of arrays stored in a structure,
%   containing fields W for weights and b for biases. It can also be a
%   neural network object.
%
% Syntax:
%    [completed] = mat2tf(model_name)
%    mat2tf(model_name)
%
% Inputs:
%    model_name: name of the .mat file which contains the NN. The input
%    might be a neural network object or a set of arrays (from NNV)
%
% Outputs:
%    completed:boolean variable which indicates if the translation
%    termianted succesfully.
%    A new NN file saved in .tf format is created automatically. It
%    contains the same name with the input model and is stored in the
%    ERAN_files/ folder.
%
% Examples:
%
% mat2tf('CartPole_Controller')
% [completed]=mat2tf('CartPolecontroller_0402')
% completed=mat2tf('CartPolecontroller_0403_tanh')
%
%
% Author:       Nikos Kekatos
% Written:      17-October-2020
% Last update:  ---
% Last revision:---

addpath(genpath('./'))

savePath=strcat('ERAN_files',filesep);

str_NN=load(model_name);

% if file name contains '.mat' ending, delete the extension.
model_name=regexprep(model_name,'.mat','');

% the models might have different structures. We look for the struct with
% the field W.
nn_object=0;
if isfield(str_NN,'net')
    if isa(str_NN.net,'network') % check class type
        net=str_NN.net;
        disp('The NN is a neural network object.')
        nn_object=1;
    end
end
if ~nn_object
    if ~isfield(str_NN,'W')
        str_NN_temp=fieldnames(str_NN);
        str_NN_temp_field=str_NN_temp{1};
        str_NN=str_NN.(str_NN_temp_field);
    end
    
    if isfield(str_NN,'W')
        W=str_NN.W;
    else
        error('The specific file has a different format and it cannot be translated!')
    end
    b=str_NN.b;
    
    if isfield(str_NN,'number_of_inputs')
        no_in=str_NN.number_of_inputs;
    else % the Weights of the first hidden layer are saved in the format [no_neurons*no_inputs]
        no_in= size(W{1},2);
    end
    
    if isfield(str_NN,'number_of_outputs')
        no_out=str_NN.number_of_outputs;
    else
        no_out= size(b{end},2);
    end
    
    if isfield(str_NN,'number_of_layers')
        no_layers=str_NN.number_of_layers;
    else
        no_layers= numel(b);
    end
    
    if isfield(str_NN,'layer_sizes')
        layer_sizes=str_NN.layer_sizes;
    else
        if  isfield(str_NN,'layer_size')
            layer_sizes=str_NN.layer_size;
        else
            for i=1:no_layers
                layer_sizes(i)=length(b{i});
            end
        end
    end
    
    %ACASXU has incorrect number of layer number. It adds the inputs
    if length(layer_sizes)~=length(b)
        clearvars layer_sizes
        for i=1:no_layers
            layer_sizes(i)=length(b{i});
        end
    end
    
    if isfield(str_NN,'activation_fcns')
        activ=cellstr(str_NN.activation_fcns);
    elseif isfield(str_NN,'act_fcns')
        activ=cellstr(str_NN.act_fcns);
    else
        % default
        for i=1:no_layers
            activ{i}='relu';
        end
    end
    % https://github.com/verivital/nnv/blob/master/code/nnv/examples/Submission/ARCH_COMP2020/benchmarks/Tora_Heterogeneous/reachTora_sigmoid.m
    if strcmp(model_name,'nn_tora_sigmoid')
        activ={'logsig','logsig','logsig','logsig'};
    elseif strcmp(model_name,'nn_tora_relu_tanh')
        activ={'relu','relu','relu','tanh'};
    end
    
    % change user_defined_activ to 1 and specify the activation
    user_defined_activ=0;
    if user_defined_activ
        warning('The user has manually added the activation functions.');
        activ={' ',' ', ' '};
    end
    
else %% net object
    if ~contains(net.name,'forward','IgnoreCase',true)
        error('The translation only works with feedforward neural networks.')
    else
        % bias
        if size(net.b,1) % 1*n cell
            b=net.b';
        else
            b=net.b;
        end
        % weights
        IW=net.IW;LW=net.LW;
        IW_cell=IW(~cellfun('isempty',IW));
        LW_cell=LW(~cellfun('isempty',LW));
        W=[IW_cell,LW_cell'];
        % no_layers
        no_layers=net.numLayers;
        % layer_sizes
        for i=1:no_layers
            layer_sizes(i)=length(b{i});
        end
        %activation functions
        for i=1:no_layers
            activ{i}=net.layers{i}.transferFcn;
        end
    end
    
end
for i=1:no_layers
    if strcmp(activ{i},'relu') || strcmp(activ{i},'poslin')
        activ{i}='ReLU';
    elseif strcmp(activ{i},'tanh') || strcmp(activ{i},'tansig')
        activ{i}='Tanh';
    elseif strcmp(activ{i},'linear') || strcmp(activ{i},'purelin')
        activ{i}='Affine';
    elseif strcmp(activ{i},'sigmoid') || strcmp(activ{i},'logsig')
        activ{i}='Sigmoid';
    end
end


fileID = fopen(strcat(savePath,model_name,'.tf'),'w');
for i_layer=1:no_layers
    %activation
    fprintf(fileID,'%s\n',activ{i_layer});
    %start weights
    fprintf(fileID,'[');
    for i_neuron=1:layer_sizes(i_layer)
        w_line=regexprep(num2str(W{i_layer}(i_neuron,:)),'\s+',',');
        %weights
        if i_neuron==layer_sizes(i_layer)
            fprintf(fileID,'[%s]',w_line);
        else
            fprintf(fileID,'[%s],',w_line);
        end
        %end weights
    end
    fprintf(fileID,']');
    % bias
    if size(b{i_layer},2)==1 % it has to be a row vector otherwise the num2str does not work
        b_line=regexprep(num2str(b{i_layer}'),'\s+',',');
    elseif size(b{i_layer},1)==1
        b_line=regexprep(num2str(b{i_layer}),'\s+',',');
    else
        error('The bias term is not correctly defined.')
    end
    fprintf(fileID,'\n[%s]\n',b_line);
    
end

fclose(fileID);
if fileID~=-1
    completed=1;
else
    completed=0;
end
end
