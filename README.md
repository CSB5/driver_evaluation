

	Author: Koh Jia Yu and Denis Bertrand
	Date: 24th Mar 2017
	Version: 1.0


What is it?


A toolbox that contains scripts to evalute cancer driver prediction methods, and compare them against prediction obtained by a collection of 18 cancer driver predictions method on 15 cancer types. 

For each cancer type, we provide ready to use genomic (SNV and CNA) and transcriptomic data.


# A database of cancer genomic and driver predictions

Evaluation data set and the associted prediction can be downloaded here: ftp://ftp2.gis.a-star.edu.sg/ConsensusDriver/evaluation_data_set.tar.gz, and uncompressed as below:

	tar -xvzf evaluation_data_set.tar.gz

The test dataset is organized as follow:

## Genomic and transcriptomic data from 15 TCGA cancer types <EVALUATION_DATA_SET/DATA>:

CNA and exome point mutation data for all cancer types were obtained from [GDAC](https://gdac.broadinstitute.org) via Firehose. 

Expression data for tumor and normal samples for all cancer types was downloaded from the [TCGA website](https://tcga-data.nci.nih.gov) (level 3). 

Samples from which the 3 data type were not present were excludued.


For each cancer type we are porviding the following files:

- GDAC_somatic_mutations.filtered.maf     
	File that contains the point mutations in maf format
- point_mutation_matrix.txt               
	Matrix (column: sample ID, row: gene name) where a cell is equal to 1 if the gene is mutated (indels, missense, nonsense and splice site variants), 0 otherwhise
- CNA_matrix.txt                          
	Matrix (column: sample ID, row: gene name) where a cell is equal to 1 if the gene have undergone a focal amplification, -1 if the gene undegone a focal deletion, 0 if the gene is consedered as copy number neutral
- normalized_expression_matrix.txt        
	Matrix (column: sample ID, row: gene name) where a cell represents the normalized expression values obtain using DEseq (see paper for detailed analysis)
- differential_expression_matrix.txt      
	Matrix (column: sample ID, row: gene name) where a cell represents the fold change obtained using DEseq by comparing each tumor to a set of normal samples (see paper for detailled analysis)

## Organization of the result files from 18 method on 15 cancer types <EVALUATION_DATA_SET/RESULT>:

The .results files for the different methods are provided for each cancer types, and are writen using the following unified format:

**Gene_name:**     HUGO gene name 

**Sample:**	       In case of methods providing sample specific prediction list of patient ID separated by ';', where the gene is predicted as driver, otherwise ALL

**Rank:**	       Rank of the gene according to the method score/p-value

**Score:**	       Score/p-value repported by the method

**Info:**	       Additional information. - if this field is empty

**Sample-specific_score:**	   In case of methods providing patient specific predictions, list of score/p-value predicted for the mutation on that patient separated by ';' otherwise keep empty


# How to evaluate your methods

1)Analyse the 15 cancer types using your method(s)

2)Convert the output file of your method in our unified result format. The file should be writen for each cancer type in EVALUATION_DATA_SET/RESULTS/CANCER_TYPE/method_name.result

3)Run our evaluation script

You may run the script <result_evaluation.pl> as below,

	perl path_to/driver_evaluation/driver_evaluation.pl --config driver_evaluation.cfg

(These directories are taken with reference to <EVALUATION_DATA_SET> directory.
The config file contains the the options listed below that are relevant to your run.

#Main directory for results

analysis_dir=full_path_to/EVALUATION_DATA_SET/RESULTS/

#Output directory for plots

final_outdir=EVALUATION_RESULT_DIRECTORY

#Script directory

script_dir=full_path_to/driver_evaluation/

#Selected method file, for additional method chosen by user          ###TO COMMENT IF PARAMETER IS NOT NEEDED

selected_method_file=extra_methods.txt

 ------------------------

--> Format of selected_method_file. Notice that the method name should be the same as the one of the one use for EVALUATION_DATA_SET/RESULTS/CANCER_TYPE/method_name.result

METHODS

<method_name 1>

<method_name 2>


eg.

METHODS

ConsensusDriver

# Output:
## Cohort level evaluation
### Concordance with gold standard
The methods were evaluated on how well their predictions identified cancer driver genes based on three standard measures: precision (fraction of predictions that belong to the gold standard), recall (fraction of the gold standard contained in the predictions) and the F1 score that combines both precision and recall. The gold standard gene lists can be found in driver_evaluation/GOLD_STANDARD/.

The following file contain the evaluation using the top 50 (top 10) predictions:

**cancer_gene_CANCER_UNION_precision_RANK_50.dat:**	 
	Matrix (column: method, row: cancer type) where a cell represents the precision

**cancer_gene_CANCER_UNION_recall_RANK_50.dat:**     
	Matrix (column: method, row: cancer type) where a cell represents the recall

**cancer_gene_CANCER_UNION_F1_RANK_50.dat:**         
	Matrix (column: method, row: cancer type) where a cell represents the F1 score

**method_name_precision_RANK_50.dat:**               
	Precision as a function of the number of predictions for method_name.

## Patient level evaluation

### Number of driver per patients
**sample_nb_driver_cat_RANK_ALL.dat:**      
	Matrix (column: method, row: number of driver [0, 1, 2-3, 4-8, 9-15, 16-25, >26) where a cell represents the fraction of patients for a number of predicted driver category

### Concordance with gold standard
The methods were evaluated on how well their patient-specific predictions identified cancer driver genes at the patient level based on three standard measures: precision, recall  and the F1 score that combines both precision and recall using the same gold standard as for the cohort based evaluation.

The following files contain the evaluation using the top 5 (top 3 and top 10) predictions:

**sample_precision_RANK_5.dat:**    
	Matrix (column: method, row: sample ID) where a cell represents the precision

**sample_recall_RANK_5.dat:**       
	Matrix (column: method, row: sample ID) where a cell represents the recall

**sample_F1_RANK_5.dat:**           
	Matrix (column: method, row: sample ID) where a cell represents the F1 score

# Prediction of actionable genes
**sample_actionable_profile_5_all.dat:**              
	Matrix (column: method, row: sample ID) where a cells represent the precision

**sample_actionable_profile_5_cancer_type.dat:**      
	Matrix (column: method, row: cancer type) where a cells represent the fraction of patients with a predicted actionable gene

## Additional files
**driver_number.dat:**        
	Matrix (column: method, row: cancer type) where a cell represents the number of drivers predicted

Contact:

If you have other questions or feedback, you may direct them to Jayce (kohjy@gis.a-star.edu.sg) and Denis (bertrandd@gis.a-star.edu.sg).



