# VADT (VCF ASE Detection Tool)


**Author**: M. Joseph Tomlinson IV  
**Last Update**: Aug 28, 2019 
**Creation Data**: May 27, 2018  

### Currently the program is in beta-phase and a manuscript drafted for submission. Please feel free to use, but please cite the github for the time period until it is published. 

## Program Description

Program takes in a raw VCF file, filters the data and then performs various statistical analysis for
detection of allele specific expression (ASE) to identify highly confident occurrences of ASE.

Overall program was originally developed in the Abasht Laboratory at the University of Delaware under the supervision of [Dr. Behnam Abasht](http://canr.udel.edu/faculty/behnam-abasht/). The statistical models and overall accuracy of their code development was overseen by [Dr. Jing Qiu](https://canr.udel.edu/faculty/jing-qiu/)

Please contact Dr. Behnam Abasht (abasht@udel.edu) with any questions or concerns. 

Important Notes:
Program was written using Python 3.6.
The input VCF file should be formatted in the following format described by GATK here:
https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it

## Program Overview

### Part 1. Filtering of Data

VADT performs various filtering steps of a VCF file before testing for ASE. The first filter performed is removing all variants that have been flagged as failing (FILTER column). Specifically the program looks for a "PASS" in the filter column to include that variant in further analysis. VADT removes all variants within a certain distance of an indel defined by the user. VADT also  removes all variants with multiple reference and alternative alleles and variants with less than the minimum quality score defined by user.

Now VADT also performs sample level filtering of the data. Sample level filtering consists of filtering based on removing all samples that are homozygous, low read count (user parameter), allele count is <1% of the total counts or "No data". So, a variant can fail for any one of these filters or can fail for a combination of these filters being implemented.

### Part 2. Detection of Reference Allele Bias

A concern with detection of ASE variants is reference allele bias where the reference genome causes the reference allele to be mapped more. So, VADT tests all "testable" variants that pass the prior filters for reference allele bias and reports the final total of (reference allele count / total allele count) for a VCF file. 

### Part 3. Binomial Test

VADT performs a binomial test on all informative "testable" biallelic samples in the dataset using SciPy's binomial test module. Utilizing the raw read counts from the VCF file. 

Link for Binomial Test: https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.binom_test.html

### Part 4 Statistical Analysis of Binomial Results

VADT performs two different types of statistical analysis of the data (varying models) to identify significant ASE results, a meta-analysis (variants) or multi-dimensional p-value adjustment (samples). Each test is examining the data in a different way either to identify significant variants or significant samples.

#### Part 4 A. Meta-Analysis (Variants)

After implementing the binomial test, VADT next performs a meta-analysis accross all tested samples on per variant basis using Fisher's Method (Fisher 1958). The underlying idea of Fisher's meta-analysis is that by combining p-values together significance may be indentified from the aggregation of samples. After getting the final p-value from the meta-analysis an FDR (Benjamini-Hochberg 1995) is implemented on the final meta p-values. The cutoff for significance of the adjusted p-values is determined by the user in the parameter file. 

Linke for Combine P-values Test: https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.combine_pvalues.html

#### Part 4 B. Multi-dimensional P-value Adjustment (Samples)

After performing the binomial test, VADT next performs a multi-dimensional p-value adjustment on a per sample basis (Guo 2010) based on a user defined p-value cutoff that gets implemented at the various stages of multi-dimensional adjustment. The multi-dimensional p-value adjustment tries to take into account the multi-diminesionality of the datasets by adjusting p-values in both the x & y direction of the data.

#### Note to Users
It is important to point out VADT is not a black box program and all intermediate files are retained for the user to help the user better understand each step performed and also identify bugs in the program. No program is perfect and can always be improved!


## Installation and Use

Download the program file and optional parameter file from website and place in local working directory on computer or server.


### Required input files:

* VCF file

Note: VCF files utilized by VADT were produced by GATK's Haplotypecaller and so VCF files produced by other programs could potentially have some bugs.

### System Requirements
**Python 3.6** and the following packages

* from __future__ import division
* import os.path
* import sys
* import time
* from scipy import stats
* from scipy.stats import combine_pvalues
* import numpy
* import copy
* from datetime import datetime
* from platform import python_version

Note: Only modules that need to be installed using "pip" or whl files are NumPy and SciPy


## Running the Program
To run the program it is possible to run in an HPC enviroment or using a parameter file. For large datasets it is highly recommended to run in an HPC enviroment because run times can be extremely long for the program  (30 minutes to 10 hours based on datasets we utilized).

Parameter File (VADT_Parameter_File.txt)

The parameter file is self-explantory and all a user needs to do is open up the file and modify the parameters defined with "--". The rest of the lines of text explaining each parameter are ignored (program properly parses). It is important the parameter file needs to be in the same directory as the program to be run.

To run the program from a HPC enivroment provide the following parameters in submission file. Feel free to change any of the below parameters based on a user's needs and data type. 

<pre>
python VADT_"version" \
--File_Name name_of_file \
--Output_File_Location location_of_desired/output \
--Indel_Exclusion_Region_Length 75 \
--Quality_Score_Minimum_for_Variants 20 \
--Minimum_Read_Counts 20 \
--Meta_BH_adj_p_value_cutoff 0.05 \
--Meta_sample_p_value_cutoff 0.05 \
--Multi_Dim_adjust_pvalue_cutoff 0.05 \
--Binomial_Probability_Value 0.5
</pre>
 
*NOTE:*
"version" refers to the current version of VADT program that was downloaded.
Slashes are important in a UNIX enviroment to tell the server "go to the next line for next parameter".

Dashed lines and the single space are important when submitting parameters. The order of the parameters does not matter, but capitalization, underscores and spelling do.

When running this program on the command line locally or on a server, beware of issues that may occur if the incorrect version of python is called or if the NumPy module is not installed properly. Some computers may automatically call a prior version of python if the command "python" is given, so a user needs to specifically call python3.6 using full pathways and also needs to install required modules (aka NumPy) to that version. NumPy installed on prior versions of python will NOT work!

*VERY IMPORTANT:* when putting pathways of the files remember to put "/" at the end of the pathway!


## Parameters 

`--File_Name` name_of_file **REQUIRED** Name of file being analyzed (give full pathway of file)

`--Output_File_Location` **REQUIRED** Location where to output the data (give full pathway of output)

Variant Level Filtering Parameters

`--Indel_Exclusion_Region_Length` (Default = 1) Length of indel exclusion region to exclude variants. Value should NOT be set to 0 because indels must be removed from datasets for ASE analysis. 

`--Quality_Score_Minimum_for_Variants` (Default = 20) Minimum quality score value for a variant

Sample Level Filtering

`--Minimum_Read_Counts` (Default = 20) Minimum total read counts (ref + alt) for a sample per variant

Meta-Analysis Cutoff Values

`-Meta_BH_adj_p_value_cutoff` 0.05 (Default = 0.05) P-value cutoff for considering an adjusted p-value significant

`--Meta_sample_p_value_cutoff` 0.05 (Default = 0.05) P-value cutoff for estimating the total tally of statistically significant samples statistically significant variants identified from the meta-analysis. 

Multi-Dimensional P-Value Adjustment (Samples)

`--Multi_Dim_adjust_pvalue_cutoff`  (Default = 0.05) P-value cutoff utilized throughout the multi-dimensional p-value adjustment algorithm.

Optional Binomial Test Probability Value

`--Binomial_Probability_Value` (Default = 0.5) Probility value for success/failure for hypothesis testing. 

Note: If a user cannot mask the reference genome, this value can be adjusted to correct for reference allele bias.  


## VADT Outputs

VADT produces a global output directory called "VADT_output_" that is date and time stamped. This is done to prevent a user from overwriting their prior result runs, which a good Computional Biologist has never done (cough cough)... Inside this global directory are three directories called: Filtering_Results, Meta_Analysis_Results and Multi_Dim_Adj_Results. Each directory will now further broken down.


### Filtering_Results

The log directory consists of 4 files that show various steps of the filtering process of VCF file.

`GATK_Failing_RNA_Seq_variants.txt` - Consists of all the variants failed by GATK (variants do not contain a PASS in the filter column). These datasets was specifically seperated from the rest of the filtering datasets because this data is huge and makes opening up the filtering file for QCing much more difficult. 

`Other_Failing_RNA_Seq_variants.txt` - Consints of all the variants failed by the various filtering criterias: located near an indel, all samples homozygous, low read count (user parameter), allele count is <1% of the total counts or "No data". The specific reason why the variant failed can be seen in the first column called "Failure" and is self explanotory.

Note: Sample_Combo_Filter refers to the fact that variant had a combination of samples fail for various reasons that ultimately resulted in the variant failing.

`Identified_Indel_Regions.txt` - A txt file that contains all the identified indels in the dataset that are used to filter the remaining data.

`Testable_Informative_Filt_Variants.txt` - A txt file that contains all passing variants that can be further analyzed. It is important to point out the original VCF format in this file has been drastically filtered/changed from the original vcf , so that each sample record now consists of four features.

Verdict: Genotype: Counts: Binomial_P_value

Verdict = Overall verdict of the sample record. Important only "Biallelic" and "Homo" tags utilized in further analysis, rest of the tags refer to why the sample failed.

Genotype: The actual genotype of the sample 0/0, 1/1, 0/1 (Only biallic variants considered)

Counts: Raw counts of the reference and alternative allele 

Binomial_P_value: P-value calculated from the binomial test of the raw counts.


### Meta_Analysis_Results (10 Files)

`Summary_Report_meta_analysis.txt` - Summary report of the entire ASE analysis

`data_for_ref_allele_bias_plotting.txt` - Counts for all testable variants, so overall reference allele bias can be investigated on a per variant basis or globally

`variant_meta_analysis_results.txt` - All meta-analysis results from the testable variants including the FDR corrected p-value and if the value is significant (user input parameter). 

`sig_meta_analysis_variants.txt` - All variants considered statistically significant from the meta-analysis

`not_sig_meta_analysis_variants.txt` - All variants not considered statistically significant from the meta-analysis

`sig_samples_report.txt` - Tallying report of samples results from the sig_meta_analysis_variants.txt file

`sig_variants_report.txt` - Tallying report of variants results from the sig_meta_analysis_variants.txt file

`maf_freq_binning_one_ase_hit.txt` - Frequency binning based minor allele frequency (maf) of all variants with atleast one ASE significant hit

`ase_freq_prevalence_among_samples.txt` - Frequency binning of ASE prevalance among all samples to see the overall occurrence of ASE among the samples

`sig_variants_report_and_meta_results.txt` - merged results from sig_variants_report.txt and variant_meta_analysis.results.txt


### Multi_Dim_Adj_Results (9 Files)

`Summary_Report_multi_dimensiona_pvalue_adj.txt` - Summary report of the entire ASE analysis

`data_for_ref_allele_bias_plotting.txt` - Counts for all the informative "testable" variants, so overall reference allele bias can be investigated on a per variant basis or globally

`multi_dim_adj_pvalues_testable.txt` - File consists of all the informative "testable" samples with a new entry in the sample column on whether that sample passed the cutoff for that specific variant. The cutoff p-value for each variant calcualted from the multi-dimensional p-value algorithm can be found in the variants format column (last entry). It is important to note each variant will have a different cutoff p-value. 

`sig_multi_dim_adj_results.txt` - All variants with atleast one sample identified as significant

`not_sig_multi_dim_adj_results.txt` - All variants that failed and did not contain one significant sample

`sig_samples_report.txt` - Tallying report of samples results from the sig_multi_dim_adj_results.txt file

`sig_variants_report.txt` - Tallying report of variants results from the sig_multi_dim_adj_results.txt file

`maf_freq_binning_one_ase_hit.txt` - Frequency binning based minor allele frequency (maf) of all variants with atleast one ASE significant hit

`ase_freq_prevalence_among_samples.txt` - Frequency binning of ASE prevalance among all samples to see the overall occurrence of ASE among the samples



### References

Benjamini, Y. and Y. Hochberg, Controlling the false discovery rate: a pratical and powerful approach to multiple testing. Journal of the Royal Statistical Society, 1995. 57(1): p. 289-300

Fisher, R., Statistical Methods for Research Workers (Thirteenth Edition-Revised). 1958, New York: Hafner Publishing Company Inc.

Guo, W., S.K. Sarkar, and S.D. Peddada, Controlling false discoveries in multidimensional directional decisions, with applications to gene expression data on ordered categories. Biometrics, 2010. 66(2): p. 485-92.








