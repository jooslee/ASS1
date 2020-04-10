# Targeting purine synthesis in ASS1 expressing cancers promotes response to immune checkpoint inhibitors
We provide the custum codes that were used to generate the results of our manuscript ‘Targeting purine synthesis in ASS1 expressing cancers promotes response to immune checkpoint inhibitors’.

## Structure of the codes
The repository includes the source code under ./R directory and relevant data under ./data directory. The code was used to perform TCGA analysis.

## Platform used to test the codes
The code was tested using R version 3.6.1 (2019-07-05) on a x86_64-pc-linux-gnu (64-bit) platform, using libraries survminer v0.4.6, survival v3.1-11, and data.table v1.12.8.

## Installing
1. Clone the git repository.

```
git clone https://github/jooslee/ASS1.git
```

2. Obtain the data from a data folder. All the relevant data is compressed in data.zip file. Once extracted it generates ./data folder, and the data for inference of the SL/SR partners and test the predictions in clinical trial datasets are available.

```
cd ASS1/data
wget https://hpc.nih.gov/~leej55/ASS1/prob.TCGA.RData
```

3. Install R libraries as needed.
```
> install.packages("data.table")
> install.packages("survminer")
> install.packages("survival")
```
4. Launch the code of interest

```
> source("./R/TCGA.R")
```


