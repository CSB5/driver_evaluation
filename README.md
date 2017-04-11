# Introduction

**Driver evaluation** is a toolbox that contains scripts to evalute cancer driver prediction methods, and compare them against prediction obtained by a collection of 18 cancer driver predictions method on 15 cancer types. 
For each cancer type, we provide ready to use genomic (SNV and CNA) and transcriptomic data.
Our database contain results from methods that differ widely in the information they require as input (e.g. point mutations, indels, CNAs, expression data etc.) and in the models/assumptions.
We classifed the methods in the folowing categegories:
- methods that belong to the functional impact category (primarily designed to identify function altering mutations but have been used for predicting drivers, such as [SIFT](http://sift.bii.a-star.edu.sg/), [PolyPhen2](http://genetics.bwh.harvard.edu/pph2/), [MutationTaster](http://www.mutationtaster.org/) and [MutationAssessor](http://mutationassessor.org/r3/)
- methods that tailor the Functional Impact idea to cancer by learning specific models such as [CHASM](http://wiki.chasmsoftware.org/index.php/CHASM_Overview), [transFIC](https://omictools.com/transformed-functional-impact-score-for-cancer-tool) and [fathmm](http://fathmm.biocompute.org.uk/) 
- methods that use cohort based analysis to detect signals of positive selection such as [ActiveDriver](http://individual.utoronto.ca/reimand/ActiveDriver/), [MutSigCV](http://software.broadinstitute.org/cancer/software/genepattern/modules/docs/MutSigCV), [OncodriveCLUST](https://bitbucket.org/bbglab/oncodriveclust) and [OncodriveFM](http://bg.upf.edu/group/projects/oncodrive-fm.php)
- methods that integrate mutation data with transcriptomic data by looking for mutation-expression correlations (OncodriveCIS)[http://bg.upf.edu/group/projects/oncodrive-cis.php] and (S2N)[http://www.maayanlab.net/S2N/]
- methods that use information from gene/protein interaction networks to analyze the effect of mutations such as [NetBox](http://cbio.mskcc.org/downloads/), [HotNet2](http://compbio.cs.brown.edu/projects/hotnet2/), [DriverNet](http://driver-net.com/about), [DawnRank](http://bioen-compbio.bioen.illinois.edu/DawnRank/) and [OncoIMPACT](https://sourceforge.net/projects/oncoimpact/). 

# A database of cancer genomic and driver predictions

Evaluation data set and the associted prediction can be downloaded here: ftp://ftp2.gis.a-star.edu.sg/ConsensusDriver/evaluation_data_set.tar.gz, and uncompressed as below:

	tar -xvzf evaluation_data_set.tar.gz

The test dataset is organized as follow:

### Genomic and transcriptomic data from 15 TCGA cancer types:

CNA and exome point mutation data for all cancer types were obtained from [GDAC](https://gdac.broadinstitute.org) via Firehose. Expression data for tumor and normal samples for all cancer types was downloaded from the [TCGA website](https://tcga-data.nci.nih.gov) (level 3). Samples from which the 3 data type were not present were excludued. These can be found in the test dataset, <EVALUATION_DATA_SET/DATA>.


For each cancer type we are porviding the following files (see paper for detailed analysis):

- **GDAC_somatic_mutations.filtered.maf**     
	File that contains the point mutations in maf format
- **point_mutation_matrix.txt**               
	Matrix (column: sample ID, row: gene name) where a cell is equal to 1 if the gene is mutated (indels, missense, nonsense and splice site variants), 0 otherwhise
- **CNA_matrix.txt**                          
	Matrix (column: sample ID, row: gene name) where a cell is equal to 1 if the gene have undergone a focal amplification, -1 if the gene undegone a focal deletion, 0 if the gene is consedered as copy number neutral
- **normalized_expression_matrix.txt**        
	Matrix (column: sample ID, row: gene name) where a cell represents the normalized expression values obtain using DEseq
- **differential_expression_matrix.txt**      
	Matrix (column: sample ID, row: gene name) where a cell represents the fold change obtained using DEseq by comparing each tumor to a set of normal samples (see paper for detailled analysis)

### Organization of the result files from 18 method on 15 cancer types:

The .result files for the different methods are provided for each cancer types. They can be found in <EVALUATION_DATA_SET/RESULT> and are writen using the following unified format:

- **Gene_name:**     HUGO gene name 

- **Sample:**	       In case of methods providing sample specific prediction list of patient ID separated by ';', where the gene is predicted as driver, otherwise ALL

- **Rank:**	       Rank of the gene according to the method score/p-value

- **Score:**	       Score/p-value repported by the method

- **Info:**	       Additional information. - if this field is empty

- **Sample-specific_score:**	   In case of methods providing patient specific predictions, list of score/p-value predicted for the mutation on that patient separated by ';' otherwise keep empty

# Dependencies

- The Perl module 'Config::Simple' is required. The steps needed to install this module can be found [here](http://www.livejournal.com/doc/server/lj.install.perl_setup.modules.html).

# How to evaluate your methods

1) Analyse the 15 cancer types using your method(s)

2) Convert the output file of your method in our unified result format. The file should be writen for each cancer type in EVALUATION_DATA_SET/RESULTS/CANCER_TYPE/method_name.result

3) Run our evaluation script

You may run the script <driver_evaluation.pl> as below,

	perl path_to/driver_evaluation.pl --config path_to/driver_evaluation.cfg

The config file contains the the options.

- **analysis_dir:**	Main directory for results

- **final_outdir:**	Output directory for plots

- **script_dir:**	Script directory

- **selected_method_file:**	Selected method file, for additional method chosen by user (To comment this if there is no additional methods)

 ------------------------

--> The file <extra_methods.txt> contains additional methods that you want to introduce to this run. 

Notice that the method name should be the same as the one of the one use for EVALUATION_DATA_SET/RESULTS/CANCER_TYPE/method_name.result


# How to perform a test run

1) The entries for the options in the config file <driver_evaluation.cfg>, in the <driver_evaluation/TEST_DATA_SET/> directory, have already been filled for your test run.
However, you may need to specify the analysis directory, analysis_dir=path_to/EVALUATION_DATA_SET/RESULTS/.

2) Run the script <driver_evaluation.pl>, in the <driver_evaluation/TEST_DATA_SET/> directory, as follows.

        cd driver_evaluation/TEST_DATA_SET/
	
	perl ../driver_evaluation.pl --config ./driver_evaluation.cfg



# Output:
## Cohort level evaluation
### Concordance with gold standard
The methods were evaluated on how well their predictions identified cancer driver genes based on three standard measures: precision (fraction of predictions that belong to the gold standard), recall (fraction of the gold standard contained in the predictions) and the F1 score that combines both precision and recall. The gold standard gene lists can be found in driver_evaluation/GOLD_STANDARD/.

The following file contain the evaluation using the top 50 predictions (similar files are also available for the top 10 predictions):

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
	Matrix (column: method, row: number of driver [0, 1, 2-3, 4-8, 9-15, 16-25, >26]) where a cell represents the fraction of patients for a number of predicted driver category

### Concordance with gold standard
The methods were evaluated on how well their patient-specific predictions identified cancer driver genes at the patient level based on three standard measures: precision, recall  and the F1 score that combines both precision and recall using the same gold standard as for the cohort based evaluation.

The following files contain the evaluation using the top 5 predictions (similar files are also available for the top 3 and top 10 predictions):

- **sample_precision_RANK_5.dat:**    
	Matrix (column: method, row: sample ID) where a cell represents the precision

- **sample_recall_RANK_5.dat:**       
	Matrix (column: method, row: sample ID) where a cell represents the recall

- **sample_F1_RANK_5.dat:**           
	Matrix (column: method, row: sample ID) where a cell represents the F1 score

### Prediction of actionable genes
- **sample_actionable_profile_5_all.dat:**              
	Matrix (column: method, row: sample ID) where a cells represent the precision

- **sample_actionable_profile_5_cancer_type.dat:**      
	Matrix (column: method, row: cancer type) where a cells represent the fraction of patients with a predicted actionable gene

## Additional files
- **driver_number.dat:**        
	Matrix (column: method, row: cancer type) where a cell represents the number of drivers predicted

# Contact:

If you have other questions or feedback, you may direct them to Jayce (kohjy@gis.a-star.edu.sg) and Denis (bertrandd@gis.a-star.edu.sg).



