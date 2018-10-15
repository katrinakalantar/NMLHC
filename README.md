# NMLHC
Noisy Metagenomically Labelled Host Classifier

This directory contains all scripts associated with the analysis and development of algorithms for host gene expression classifiers in the context of noisy metagenomic labels.



### Classification with Noisy Labels: /robust_logisticregression 
Contains matlab implementation of the rLR algorithm published by [Bootkrajang _et al_](http://www.cs.science.cmu.ac.th/person/jakramate/) in 2012 for inferring logistic regression classifiers from data with noisy labels.

The code is run through R.matlab, which opens a server to enable execution of matlab code from R scripts. This is a first-pass implementation to benchmark the methods rapidly.


### Parameters 

The scripst are written to take in a **parameters.json** file containing all information required for running that experiment.

Here is an example of the overall .json structure:

```

{
  "array": [
    1,
    2,
    3
  ],
  "ITER": 10,
  "CLS": 2,
  "DATASET_PARAM": "generate_data",
  "use_PCs": false,
  "feature_select": true,
  "common_reg": "lasso",
  "common_sn": 1e-8,
  "common_maxIter": 1000,
  "geneset_file": "/Users/kkalantar/Documents/Research/NMLHC/reference/c2.cp.v6.2.symbols.gmt",
  "convert_genenames": true,
  "collapse_pathways": false,
  "datasets":[
    {
      "name": "BACTERIA_VIRUS",
      "type": "geo",
      "series_filename": "/Users/kkalantar/Documents/Research/NMLHC/Exp1_HostBenchmark/data/GSE60244_series_matrix.txt",
      "source_variable": "source_name_ch1"
    }
  ],
  "DS_SIZE_RANGE":[
    200
  ],
  "DIM_RANGE":[
    2000
  ],
  "EXP_RANGE":[
    0.0,
    0.1,
    0.2,
    0.3
  ],
  "EXP_RANGE_J":[
    0.0,
    0.2,
    0.3
  ],
  "estG": false
}

```


Here are dataset-specific .json blocks for pre-curated datasets used in the analysis:

[**GSE60244**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60244) microarray dataset from _Transcriptional Profiling is Superior to Procalcitonin to Discriminate Bacterial vs. Viral Lower Respiratory Tract Infections in Hospitalized Adults_ where whole blood was collected from individuals with Bacterial, Viral, or coinfection.

```
{
      "name": "BACTERIA_VIRUS",
      "type": "geo",
      "series_filename": "/Users/kkalantar/Documents/Research/NMLHC/Exp1_HostBenchmark/data/GSE60244_series_matrix.txt",
      "source_variable": "source_name_ch1"
}
```

[**GSE81046**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81046) RNA-seq dataset from _Genetic ancestry and natural selection drive population differences in immune responses to pathogens in humans_ where primary macrophages were isolated from individuals of differing ancestry and exposed to live bacterial pathogens in culture.

```
{
      "name": "Listeria_Non-infected",  # can also use Salmonella as a contrast here
      "type": "geo",
      "series_filename": "/Users/kkalantar/Documents/Research/NMLHC/geo_data/ANCESTRY_DATASET.rds",
      "source_variable": "infection:ch1"
}
```

[**GSE92904**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92904) RNA-seq dataset from _Genetic analysis of isoform usage in the human anti-viral response_ in which monodyte-derived dendritic cells were isolated from individuals and exposed to viral or IFN-beta treatment.

```
{
      "name": "IFN-beta_baseline",  # can also use influenza as a contrast here
      "type": "geo",
      "series_filename": "/Users/kkalantar/Documents/Research/NMLHC/geo_data/ANTIVIRAL_DATASET.rds",
      "source_variable": "characteristics_ch1.3"
}
```


