### Developing Structure-Activity Relationship (SAR): 

The primary objective of QSAR is to establish a relationship between the chemical structure of compounds and their corresponding activities. This relationship can be developed using both traditional ML algorithms and more advanced DL techniques.

#### ML-based QSAR models :
- typically involve feature selection or extraction, followed by training and evaluation of ML algorithms such as linear regression, decision trees, random forests, support vector machines, or gradient boosting. These models aim to identify the most relevant molecular descriptors or fingerprints that contribute to the observed activity and use them to make predictions for new compounds.
- Here, I tried exploring 3 different features :
    
    a. ECFP / Morgan finger prints
    
    b. MACCS
    
    c. Physiochemical properties
    
    with Random forest regressor to understand how different features plays a good role in the model building.
    
 - Also, I have explored using Pycaret AutoML.
    
#### DL-based QSAR models :
-leverage the power of deep neural networks to learn complex representations directly from the molecular structures. Convolutional neural networks (CNNs) or graph convolutional networks (GCNs) can be employed to capture the hierarchical features of molecules. DL models can handle raw molecular structures and learn both local and global structural patterns, potentially leading to more accurate predictions.
- Here, I have tried exploring:

A. Simple FNN using MACCS fingerprint
B. RNN-LSTM using mol2vec 
- Thing will be exploring soon:

1. MPNN / Message passing neural network
2. GAN / Graph attention network
3. GCN / Graph convolutional networks
4. AFP / Attentive fingerprint (AFP)
