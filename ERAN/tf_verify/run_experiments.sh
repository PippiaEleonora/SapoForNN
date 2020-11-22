
printf "\n--------------------------------------------------------------\n"
printf "\n\t\t Starting the experiments\n"
printf "\n--------------------------------------------------------------\n"

read -p "Press <ENTER> to start Experiments:"
printf "\n"
read -p "Choose model, 1 (synthetic), 2 (NNV toy example) or 3 (NNV ACC model): "
if [[ $REPLY =~ ^[1]$ ]];then
    model="exp555.tf"
    no=1
    args=" "
else 
    if [[ $REPLY =~ ^[2]$ ]];then
        model="Engine_Toy_Tansig_net.tf"
        no=2
        args="--num_inputs 2"
    fi
fi
if [[ $REPLY =~ ^[3]$ ]];then
    model="ACC_controller_3_20_tanh.tf"
    no=3
    args="--num_inputs 5"
fi

printf "\nThe selected model is %s. \n" "$model"
printf "The number ID is %s.\n\n" "${no}"
if [[ $no = 1 ]];then
    printf "\n \tRunning default/original experiment"
    printf "\n Default configuration"
    printf "\n model: exp555.tf"
    printf "\n input: [-1,1]"
    printf "\n sparse_n: 3"
    printf "\n splitting: yes"
    printf "\n sanity_check: yes"
    printf "\n polynomial: hard coded\n\n"
else 
if [[ $no = 2 ]]; then
printf "\n \tRunning default/original experiment"
printf "\n Default configuration"
printf "\n model: Engine_Toy_Tansig.tf"
printf "\n input: [-1,1] x [-1.1]"
printf "\n sparse_n: 3"
printf "\n splitting: yes"
printf "\n sanity_check: yes"
printf "\n polynomial: hard coded\n\n"
else
if [[ $no = 3 ]]; then
printf "\n \tRunning default/original experiment"
printf "\n Default configuration"
printf "\n model: ACC_controller_3_20_tansig.tf"
printf "\n input: [-1,1]^5"
printf "\n sparse_n: 3"
printf "\n splitting: yes"
printf "\n sanity_check: yes"
printf "\n polynomial: hard coded\n\n"
fi
fi
fi
read -p "Do you want to continue with experiment 1A (y or n)?"
if [[ $REPLY =~ ^[Yy]$ ]]
then
python3 . --netname DATA/$model --domain refinepoly --epsilon 0.1 --sparse_n 2 \
--refine_neurons --num_test 1 $args
printf "\n\n"
fi
read -p "Do you want to continue with experiment 1B (y or n)?"
if [[ $REPLY =~ ^[Yy]$ ]]
then
printf "\n\t Experiment 1B:  evaluate the impact of sanity check.\n"
printf "\nflag --sanity_check 0"
printf "\nrest: default\n\n "
sleep 2
python3 . --netname DATA/$model --domain refinepoly --epsilon 0.1 --sparse_n 2 \
--refine_neurons --num_test 1 --sanity_check 0 $args
printf "\n\n"
printf ""
fi
read -p "Do you want to continue with experiment 2 (y or n)?"
if [[ $REPLY =~ ^[Yy]$ ]]
then
printf "\n\t Experiment 2A:  different input ranges\n"
printf "\ninput: [-2,2]"
printf "\nrest: default\n\n "
sleep 2
python3 . --netname DATA/$model --domain refinepoly --epsilon 0.1 --sparse_n 2 \
--refine_neurons --num_test 1 --lower_bound -2 \
--upper_bound 2 $args
printf "\n\n"
fi
read -p "Do you want to continue with experiment 2B (y or n)?"
if [[ $REPLY =~ ^[Yy]$ ]]
then
printf "\n\t Experiment 2B:  different input ranges"
printf "\ninput: [-5,5]"
printf "\nrest: default\n\n"
sleep 2
python3 . --netname DATA/$model --domain refinepoly --epsilon 0.1 --sparse_n 2 \
--refine_neurons --num_test 1 --lower_bound -5 \
--upper_bound 5 $args
printf "\n\n"
fi

read -p "Do you want to continue with experiment 2C (y or n)?"

if [[ $REPLY =~ ^[Yy]$ ]]
then
printf "\t Experiment 2C:  different input ranges"
printf "\ninput: [-10,10]"
printf "\nrest: default\n\n"
sleep 2
python3 . --netname DATA/$model --domain refinepoly --epsilon 0.1 --sparse_n 2 \
--refine_neurons --num_test 1 --lower_bound -10 \
--upper_bound 10  $args
printf "\n\n"
fi

read -p "Do you want to continue with experiment 3A (y or n)?"
if [[ $REPLY =~ ^[Yy]$ ]]
then
printf "\n\t Experiment 3A:  different grouping term k-tanh\n"
printf "\nsparse_n: 1"
printf "\ninputs: [-5,5]"
printf "\nrest: default\n\n "
sleep 2
python3 . --netname DATA/$model --domain refinepoly --epsilon 0.1 --sparse_n 1 \
--refine_neurons --num_test 1 --sanity_check 1 --splitting 1 --lower_bound -5 \
--upper_bound 5 $args
fi

read -p "Do you want to continue with experiment 3B (y or n)?"
if [[ $REPLY =~ ^[Yy]$ ]]
then
printf "\n\t Experiment 3B:  different grouping term k-tanh\n"
printf "\nsparse_n: 2"
printf "\ninputs: [-5,5]"
printf "\nrest: default\n\n "
sleep 2
python3 . --netname DATA/$model --domain refinepoly --epsilon 0.1 --sparse_n 2 \
--refine_neurons --num_test 1 --sanity_check 1 --splitting 1 --lower_bound -5 \
--upper_bound 5 $args
fi
read -p "Do you want to continue with experiment 3C (y or n)?"
if [[ $REPLY =~ ^[Yy]$ ]]
then
printf "\n\t Experiment 3C:  different grouping term k-tanh\n"
printf "\nsparse_n: 3"
printf "\ninputs: [-5,5]"
printf "\nrest: default\n\n "
sleep 2
python3 . --netname DATA/$model --domain refinepoly --epsilon 0.1 --sparse_n 3 \
--refine_neurons --num_test 1 --sanity_check 1 --splitting 1 --lower_bound -5 \
--upper_bound 5 $args
fi

read -p "Do you want to continue with experiment 4 (y or n)?"
if [[ $REPLY =~ ^[Yy]$ ]]
then
printf "\n\t Experiment 4A:  splitting vs no-splitting\n"
printf "\nrest: default\n\n "
printf "\nSplitting: True or False\n\n"
sleep 2
python3 . --netname DATA/$model --domain refinepoly --epsilon 0.1 --sparse_n 3 \
--refine_neurons --num_test 1 --splitting 0 --lower_bound -1 \
--upper_bound 1 $args
printf "\n The results without splitting are displayed above. The same experiment was done"
printf "\n with splitting before (Experiment 1)."
if [[ $no = 1 ]]; then
printf  "The results are shown below for comparison. "
printf "\n The total time ERAN is  0.004237651824951172 seconds"
printf "\n lower bounds _eran:[-4.448640884943766]"
printf "\n upper bounds _eran:[6.90203888029184]"
printf "\n The total time ERAN+SAPO is  144.8912649154663 seconds"
printf "\n lower bounds _sapo:[-4.37780926261568]"
printf "\n upper bounds _sapo:[6.876765039131466]\n\n"
fi
fi 

read -p "Do you want to continue with experiment 5 (y or n)?"
if [[ $REPLY =~ ^[Yy]$ ]]
then
printf "\n\t Experiment 5A:  polynomial approx.\n"
printf "\nrest: default\n\n "
printf "\nSplitting: True or False\n\n"
sleep 2
#python3 . --netname DATA/$model --domain refinepoly --epsilon 0.1 --sparse_n 3 \
#--refine_neurons --num_test 1 --splitting 1 --sanity_check 1 --lower_bound -1 \
#--upper_bound 1 --poly_dynamic 1 --poly_order 5
fi
printf "\n Experiment 5 uses the polynomial approximation created in Python and called in Sapo."
printf "\n The error is Gdk-Message: 21:17:28.669: .: Fatal IO error 0 (Success) on X server :0."
printf "\n\n This is the end of the experients. \n"
