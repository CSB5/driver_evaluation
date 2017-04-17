# Introduction

This repository contains scripts to evaluate cancer driver prediction methods and compare them against predictions obtained from 18 different existing methods on 15 TCGA cancer types. For each cancer type, we provide ready to use genomic (point mutations and copy number variation) and transcriptomic data. In addition, we provide a database containing results from the different driver prediction methods ([SIFT](http://sift.bii.a-star.edu.sg/), [PolyPhen2](http://genetics.bwh.harvard.edu/pph2/), [MutationTaster](http://www.mutationtaster.org/), [MutationAssessor](http://mutationassessor.org/r3/), [CHASM](http://wiki.chasmsoftware.org/index.php/CHASM_Overview), [transFIC](https://omictools.com/transformed-functional-impact-score-for-cancer-tool), [fathmm](http://fathmm.biocompute.org.uk/), [ActiveDriver](http://individual.utoronto.ca/reimand/ActiveDriver/), [MutSigCV](http://software.broadinstitute.org/cancer/software/genepattern/modules/docs/MutSigCV), [OncodriveCLUST](https://bitbucket.org/bbglab/oncodriveclust), [OncodriveFM](http://bg.upf.edu/group/projects/oncodrive-fm.php), [OncodriveCIS](http://bg.upf.edu/group/projects/oncodrive-cis.php), [S2N](http://www.maayanlab.net/S2N/)[NetBox](http://cbio.mskcc.org/downloads/), [HotNet2](http://compbio.cs.brown.edu/projects/hotnet2/), [DriverNet](http://driver-net.com/about), [DawnRank](http://bioen-compbio.bioen.illinois.edu/DawnRank/) and [OncoIMPACT](https://sourceforge.net/projects/oncoimpact/)). For further details see [Bertrand et al, 2017]. 

# Database containing cancer genomic data and driver predictions

Evaluation datasets and predictions for the 18 different methods listed above. The database will be automatically downloaded during the software installation (see below).
The evaluation dataset is organized as follows:

### Genomic and transcriptomic data from 15 TCGA cancer types

Point mutation and copy number variation data for all cancer types was obtained from [GDAC](https://gdac.broadinstitute.org) via Firehose. Expression data for tumor and normal samples for all cancer types was downloaded from the [TCGA website](https://tcga-data.nci.nih.gov) (level 3). Samples for which the 3 data types were not available were excluded. 

The directory `EVALUATION_DATA_SET/DATA` contains the following files for each of the cancer types:

- **GDAC_somatic_mutations.filtered.maf**     
	File that contains point mutations in maf format
- **point_mutation_matrix.txt**               
	Matrix (column: sample ID, row: gene name) where a cell is equal to 1 if the gene is mutated (indels, missense, nonsense and splice site variants), 0 otherwise
- **CNA_matrix.txt**                          
	Matrix (column: sample ID, row: gene name) where a cell is equal to 1 if the gene is part of a focal amplification, -1 if the gene is part of a focal deletion, 0 otherwise
- **normalized_expression_matrix.txt**        
	Matrix (column: sample ID, row: gene name) where a cell represents the normalized expression value obtained using DESeq
- **differential_expression_matrix.txt**      
	Matrix (column: sample ID, row: gene name) where a cell represents the fold change obtained using DESeq by comparing each tumor to a set of normal samples (see [Bertrand et al, 2017])

### Predictions from 18 methods on 15 cancer types

The .result files for the different methods are provided for each cancer type. They can be found in `EVALUATION_DATA_SET/RESULT` and are in the following unified format:

- **Gene_name:**     HUGO gene name 

- **Sample:**	       In the case of methods that provide patient specific predictions, list of patient IDs where the gene is predicted as a driver (separated by ';'), otherwise ALL

- **Rank:**	       Rank of the gene according to the method based on the reported score or p-value

- **Score:**	       Score or p-value reported by the method

- **Info:**	       Additional information reported by the method

- **Sample-specific_score:**	   In the case of methods that provide patient specific predictions, list of scores or p-values provided by the method per patient (separated by ';')

# Installing the evaluation scripts

- Download the latest version of the software and unzip it or `git clone https://github.com/CSB5/driver_evaluation.git`

- Run the `./install.pl` command to install and download the required databases.

# Evaluate new methods using the evaluation scripts

1) Analyze data for the different cancer types using the new method

2) Convert the output file of the method into the unified result format defined above. The result files should be named according to the cancer type analysed with a '.result' extention (e.g. GBM.result) and organized in a single directory. The name of the directory will be used as the method name in the evaluation result files.

3) Run the evaluation script `driver_evaluation.pl` using the following options:

   - **--method_dir:**	Directory that contains the '.result' files of the evaluated method

   - **--out_dir:**	The directory that will contain the evaluation result files

You can perfom a test run using the following commands:

~~~~
cd TEST_DATA_SET/
perl ../bin/driver_evaluation.pl --method_dir ConsensusDriver/ --out_dir eval_result
~~~~

------------------------


# Description of the evaluation result files
## Cohort level evaluation
### Concordance with gold standard
The methods were evaluated on how well their predictions identified cancer driver genes based on three standard measures: precision (fraction of predictions that belong to the gold standard), recall (fraction of the gold standard contained in the predictions) and the F1 score that combines both precision and recall. The gold standard gene lists can be found in `driver_evaluation/GOLD_STANDARD/`.

The following files contain the evaluation results using the top 50 predictions (similar files are also available for the top 10 predictions):

- **cancer_gene_CANCER_UNION_precision_RANK_50.dat:**
	Matrix (column: method, row: cancer type) where a cell represents the precision

- **cancer_gene_CANCER_UNION_recall_RANK_50.dat:**     
	Matrix (column: method, row: cancer type) where a cell represents the recall

- **cancer_gene_CANCER_UNION_F1_RANK_50.dat:**         
	Matrix (column: method, row: cancer type) where a cell represents the F1 score

- **method_name_precision_RANK_50.dat:**               
	Precision as a function of the number of predictions for method_name.

## Patient level evaluation

### Number of driver per patients
- **sample_nb_driver_cat_RANK_ALL.dat:**      
	Matrix (column: method, row: number of predicted driver category [0, 1, 2-3, 4-8, 9-15, 16-25, >26]) where a cell represents the fraction of patients for which the predicted number of drivers falls in the given category.

### Concordance with gold standard
The methods were evaluated on how well their patient-specific predictions identified cancer driver genes at the patient level based on three standard measures: precision, recall and the F1 score.

The following files contain the evaluation results using the top 5 predictions (similar files are also available for the top 3 and top 10 predictions):

- **sample_precision_RANK_5.dat:**    
	Matrix (column: method, row: sample ID) where a cell represents the precision

- **sample_recall_RANK_5.dat:**       
	Matrix (column: method, row: sample ID) where a cell represents the recall

- **sample_F1_RANK_5.dat:**           
	Matrix (column: method, row: sample ID) where a cell represents the F1 score

### Prediction of actionable genes
The method were evaluated on their ability to identify patient specific drivers that are potentially actionable. The actionable gene list was obtained by combining actionable gene lists from [intOGen](https://www.intogen.org/downloads) and [OncoKB](http://oncokb.org/#/), and can be found in `driver_evaluation/ACTIONABLE_GENES/combine_target.dat`. The following files contain the evaluation results using the top 5 predictions (similar files are also available for top 10 predictions):

- **sample_actionable_profile_5_all.dat:**              
	Matrix (column: method, row: actionable gene category [0: approved drug, 1: investigational target, 2: research target, 3: not actionable]) for each patient we retain the predicted driver that falls in the best actionable gene category (with 0 being the best one and 3 the worst), and each cell represents the fraction of patients that falls in a given actionable gene category.

- **sample_actionable_profile_5_cancer_type.dat:**      
	Matrix (column: method, row: cancer type) where a cell represents the fraction of patients with a predicted actionable gene.

## Additional files
- **driver_number.dat:**        
	Matrix (column: method, row: cancer type) where a cell represents the number of drivers predicted

# Contact:

Please direct any questions or feedback to Denis Bertrand (bertrandd@gis.a-star.edu.sg) and Niranjan Nagarajan (nagarajann@gis.a-star.edu.sg).
