# SapoForNN

## Download
Open your terminal, navigate to your desired directory and write:

``git clone https://github.com/PippiaEleonora/SapoForNN.git``

## Installation

To run the following project you need to install all the dependences of the [Sapo](https://github.com/dreossi/sapo) project and you need to install the [ERAN](https://github.com/eth-sri/eran/tree/45edbce4dcbeeffb1d77c4f056f2725868b73ef5) tool separately.

Check that your ERAN folder has both the ELINA folder and the "data" folder inside.
You need to define the ERAN path
``cd ERAN/tf_verify``

Open the '\_\_main\__.py' and define the variable 'Eranfolder', at line 4, with your already installed ERAN (wher you have ELINA, deepG etc.)


To compile the Sapo project:
``cd Sapo
cmake -DBUILD_SHARED_LIBS=on
make``

## Example
Run the first example
`` python3 . --netname DATA/mnist_relu_3_50.tf --epsilon 0.1 --domain refinepoly --dataset mnist --sparse_n 3 --refine_neurons --num_test 1``

