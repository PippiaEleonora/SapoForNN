# SapoForNN

## Download
Open your terminal, navigate to your desired directory and write:
``git clone https://github.com/PippiaEleonora/SapoForNN.git``

## Installation

To run the following project you need to install all the dependences of the [Sapo](https://github.com/dreossi/sapo) project and you need to install the [ERAN](https://github.com/eth-sri/eran/tree/45edbce4dcbeeffb1d77c4f056f2725868b73ef5) tool.

Check that your ERAN folder has both the ELINA folder and the "data" folder inside.
You need to define the ERAN path.

Open 'eranpath.txt' that you find inside 'SapoForNN/ERAN/tf_verify'. Define the path of your already installed ERAN (where you have ELINA, deepG etc.)


To compile the Sapo project:

``` 
cd Sapo 
cmake -DBUILD_SHARED_LIBS=on
make
```


## Example
Run the first example

``` 
cd ERAN/tf_verify
python3 . --netname DATA/mnist_tanh_3_5.tf --epsilon 0.1 --domain refinepoly --sparse_n 3 --refine_neurons
```

