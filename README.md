# Modeling-Microbial-Interactions
Code from my Thesis! 

## Project Description

The project was developed as part of my thesis in collaboration with Dr. Gregory Matthews, Dr. Catherine Putonti, and Dr. George Thiruvathukal.

Bioinformaticians in the microbiology field are often interested in analyzing the microbial profiles of samples in diseased or infected cohorts. Taxonomic profiling through metagenomic sequencing can enable the recognition of such microbes or interactions between microbes that could increase risk for some outcome. For example, it may be the presence of microbes A AND B, AND NOT microbe C that strongly increase the odds of disease. It may also be the presence of microbes C OR D that strongly decrease the odds of the disease. If we could assign a weight to each of these Boolean conditions, which encapsulate interactions between different microbes, we could identify the main contributors to specific diseases. Logic regression is a statistical model that is able to quantitatively tie such Boolean statements back to changing the odds of an outcome variable. However, logic regression only takes in binary variables as features — thus the continuous relative abundances must be converted to binary variables. In order to binarize the relative abundances, we can generate a threshold vector (1xn), where n is the number of features in the data. Any element (i.e., microbial abundance) below its respective element in the threshold vector is inputted as 0, and any abundance above its respective element in the threshold vector is inputted as 1, as shown below. 
          <img width="626" height="168" alt="image" src="https://github.com/user-attachments/assets/df9765e6-c2ff-4e9e-9179-48c4ae6743c9" />

However, this threshold vector is unknown: we do not know the abundance of taxa needed to be considered present. To find this threshold vector, we use simulated annealing to tune the threshold vector. At each iteration, the threshold vector is slightly tuned (changed), the data is binarized with that particular threshold vector, and then this binarized data is used to fit a logic regression model to estimate the loss using cross validation. After finding a threshold vector that minimizes the loss through simulated annealing, the
threshold vector is used to binarize the entire data set, a logic regression model is fit, and the loss and accuracy is measured. This threshold vector also has biological relevance – it tells us the relative abundance needed in order for a microbe to be present in the model, and therefore its ability to contribute to an increase/decrease in the odds of the outcome variable.

Since many diseases often manifest as an imbalance of commensal and pathogenic taxa in the microbiota (dysbiosis), it is not enough to know who is present. To complete the picture, we must consider who is dominating the microbiota, thus necessitating an analysis of the relative abundances. Identifying abundance thresholds can improve diagnostic modeling by quantifying the levels of specific taxa needed to associate with infection or disease onset. 

## Getting Started

### Dependencies
* dplyr
* tidyverse
* caret
* LogicReg
* Metrics
* foreach
* parallel
* combinat


## Authors
Niru Shanbhag (shanbhagn123@gmail.com)


## Acknowledgments
* [Logic Reg Paper] (https://research.fredhutch.org/kooperberg/en/software/logic.html)

