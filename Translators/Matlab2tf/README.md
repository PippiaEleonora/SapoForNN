Neural Network Translator from MATLAB to ERAN format
===

In this folder, we provide a script to translate Neural Networks saved as MATLAB structures into ERAN *tf* format. We have also added support for the translation of network objects (network structure) to ERAN tf format. 

Files
---
There is a script to load and translate the Matlab NNs. The folder *MATLAB_files* contains the NN models taken from NNV [repo](https://github.com/verivital/nnv). The folder *ERAN_files* contained the translated models.

Usage
---
In MATLAB, navigate to the `Translators/matlab2tf` folder run in the command window the function
```mat2tf(model_name).``` Note that you can skip adding the file type, i.e. .MAT. For example,  you can run ```mat2tf('CartPole_Controller')``` and ``mat2tf('CartPole_Controller.mat')``.

> You can run the command `>> translate_all` to translate **all** the NN models in the *MATLAB_files* folder. 

Remarks - ERAN syntax
---
ERAN has its own format to define NNs which they call .tf. It relies on the numpy syntax of arrays complemented by the name of the activation function (supported functions: ReLU, Tanh, Sigmoid, Affine). It goes layer-by-layer and then neuron-by-neuron. Ignoring normalisation parameters and focusing on feedforward NNs, the order is as follows, 1st line activation function, 2nd line weights, 3rd line biases. 

Models - Missing Info
---

Most of the original models (NNs) are taken from the [NNV repository](https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission).
There are models that they do not have their activation functions specified. 

- The `nn_tora_relu_tanh` contains 3 hidden layers with ReLU and the output layer has a tanh. Info [here](https://github.com/verivital/nnv/blob/master/code/nnv/examples/Submission/ARCH_COMP2020/benchmarks/Tora_Heterogeneous/Specifications.txt).
- All the activation functions for the `nn_tora_sigmoid` are sigmoid. Info [here](https://github.com/verivital/nnv/blob/master/code/nnv/examples/Submission/ARCH_COMP2020/benchmarks/Tora_Heterogeneous/Specifications.txt).
- The user can manually specify the activation functions via the `user_defined_activ` flag (lines 119-123) 

MATLAB Models - Not only ReLU 
---

>[NNV](https://github.com/verivital/nnv)


1. `nn_tora_relu_tanh`

2. `nn_tora_sigmoid`

3. `ACC_controller_3_20_tanh`

4. `CartPolecontroller_0403_tanh`

5. `Engine_Toy_Tansig_net` Note that this model is not currently supported as it employs a `net` object.

6. `ACASXU_run2a_2_4_batch_2000`
7.  `mnist5x50`
8. `good_controller` neural network object (from our work)
9. `Engine_Toy_Tansig_net` neural network object (NNV)

ERAN Models -- Not only ReLU
---

>[ERAN](https://github.com/eth-sri/eran)

1. `ffnnSIGMOID__Point_6_500.pyt`

2. `ffnnTANH_Point_6_500.pyt`

Validating the translation
---
- Run the command 

```
cd SapoForNN/tf_verify
python3 . --netname mnist5x50.tf --dataset mnist --domain deepzono
```
- Run different examples without images
  
To-DO
---

1. Change the `__main__.py` in order to perform regression insted of classification.