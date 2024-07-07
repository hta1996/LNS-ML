# LNS-ML
Anytime machine learning-guided large neighborhood search for multi-agent path finding.

This code is built upon the code for MAPF-LNS.


## Setup

Follow instructions in [MAPF-LNS](https://github.com/Jiaoyang-Li/MAPF-LNS) to compile the c++ code in the folder MAPF-ML-LNS/.


Download the MAPF instances from the [MAPF benchmark](https://movingai.com/benchmarks/mapf/index.html) and place all the instances and maps in the instances/ folder. In the folder, all maps should be placed in subfolder mapf-map/ and all instances should be placed in subfolder scen-random/. 


## Training

Run the training script in the training folder to train a model for the random-32-32-10 map with 200 agents:

```
python training/train.py
```

You could provide random seed and number of iterations using options ```--seed``` and ```--iter```.


## Testing


Once finished training, there will be a trained SVM rank model in ```training/saved_model```. Open the file and extract the weights into a new weight file ```weight_file``` with the content being rows of pairs of ```featureID weights```.

Run testing with the example command:

```
./lns -m instances/mapf-map/random-32-32-10.map -a instances/scen-random/random-32-32-10-random-1.scen -k 200 -t 60 --rand=0 --lbns=5 --ubns=16 --maxIterations=1000000 --initAlgo=PPS --destoryStrategy=Adaptive --replanAlgo=PP -o test_random-32-32-10-200.csv --stats=test_random-32-32-10-200_stats.csv --samples=20 --numModels=1 --testing 1 --weight weight_file
```



# Citation

Please cite our paper if you use this code in your work.
```
@inproceedings{huang2022anytime,
  author    = {Taoan Huang and Jiaoyang Li and Sven Koenig and Bistra Dilkina},
  title     = {Anytime Multi-Agent Path Finding via Machine Learning-Guided Large Neighborhood Search},
  booktitle = {Proceedings of the AAAI Conference on Artificial Intelligence (AAAI)},
  pages     = {9368--9376},
  year      = {2022}
}
```
