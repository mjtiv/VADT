#!/usr/bin/env python3.6

"""

PROGRAM: VCF ASE Detection Tool (VADT)
DESCRIPTION: Program takes in a raw VCF file, filters the data and then performs analysis for
allele specific expression (ASE) with various statistical corrections to identify highly confident
ASE variants and samples.

AUTHOR: M. Joseph Tomlinson IV

Program Originally Created: 5-23-2018

Last Updated: 2-05-2021

"""

# Required Modules
from __future__ import division
import os.path
import sys
import time
from scipy import stats
from scipy.stats import combine_pvalues
import numpy as np
import copy
from datetime import datetime
from platform import python_version


#########################################################################################################################
############################################ Help Menu ##################################################################
#########################################################################################################################

#Help Menu for code
def help_menu_prompt ():
    print ("Help Menu\n")
    print ("To Run Program Please Supply the Following Required Parameters\n")
    print ("")
    print ("Required Parameter")
    print ("--File_Name Name_of_File.vcf")
    print ("Note: Give full pathway otherwise program assumes local directory")
    print ("")
    print ("")
    print ("Optional Parameters")
    print ("")
    print ("--Output_File_Location (default setting output to local directory)")
    print ("")
    print ("Variant Level Filtering")
    print ("--Indel_Exclusion_Region_Length (default setting = 0)")
    print ("--Quality_Score_Minimum_for_Variants (default setting = 20)")
    print ("")
    print ("Sample Level Filtering")
    print ("--Minimum_Read_Counts (default setting = 20)")
    print ("")
    print ("Meta-Analysis Cutoff Values")
    print ("--Meta_BH_adj_p_value_cutoff (default setting = 0.05)")
    print ("--Meta_sample_p_value_cutoff (default setting = 0.05) ")
    print ("")
    print ("Multi-Dimensional P-Value Cutoff Value")
    print ("--Multi_Dim_adjust_pvalue_cutoff (default setting = 0.05)")
    print ("")
    print ("Binomial Probability Value")
    print ("--Binomial_Probability_Value (default setting = 0.5)")
    print ("")
    print ("")
    print ("")
    print ("The order of the parameters does not matter, but capitalization, underscores and spelling does matter")
    print ("VERY IMPORTANT when putting pathways of the file \\ at the end is extremely important!!!")
    print ("")
    print ("Definition of Parameters")
    print ("File_Name - Name of the file being analyzed, give full pathway otherwise program assumes local directory")
    print ("Output_File_Location - Location to write data ---REMEMBER forward slash /")
    print ("Indel_Exclusion_Region_Length - Length of region around indels (forward & reverse) to filter variants.")
    print ("Quality_Score_Minimum_for_Variants - Variant quality score (phred score) minimum to filter variants ")
    print ("")
    print ("Meta-Analysis Cutoff Values")
    print ("Meta_BH_adj_p_value_cutoff - Meta-Analysis BH adjusted p-value cutoff for significance")
    print ("Meta_sample_p_value_cutoff - P-value Used for tallying samples after identification of significant variants.")
    print ("     IMPORTANT: This technique is only for a estimation of sample counts, the multi-dimensional p-value adjustment")
    print ("                gives a more accurate value")
    print ("")
    print ("Multi-Dimensional P-Value Adjustment")
    print ("Multi_Dim_adjust_pvalue_cutoff - P-value threshold for calculating significance")
    print ("")
    print ("Binomial Test Probability Value")
    print ("Binomial_Probability_Value - Probability value for testing statistical deviation.  /")
    print ("    Used for correcting reference allele bias.")
    print ("")
    print ("EXAMPLE SUBMISSION IN UNIX")
    print ("")
    print ("python name_of_program \\")
    print ("--File_Name Blue_Chickens.vcf \\")
    print ("--Indel_Exclusion_Region_Length 75 \\")
    print ("--Output_File_Location your/favorite/output/directory/ \\")
    print ("")
    sys.exit()


#########################################################################################################################
##################################### Code to Parse the Parameter File ##################################################
#########################################################################################################################


def test_Number_Input(value):

    '''

    Tests input from the user to see if its a number (decimals acceptable).
    Converts any floats to numbers and tests using the python built-in types "isdgit"
    Note: isdigit does NOT accept decimal places, so input needs to be slightly modified)
    
    : Param value: Value being tested

    : Return result: Pass or fail if its a number or not

    '''

    # Replace all periods with nothing and test if numbers
    test=value.replace('.','').isdigit()
    
    # Examines the results (True/False)
    if test==True:
        result='Pass'
    else:
        result='Fail'
 
    #return the results for the function
    return (result)


def parsing_input_parameter_file(program_parameters):

    """

    Parses apart the VADT_Parameter_File.txt to get all user input parameters
    for analyzing the data. Function also prints all user parameters to
    command line, so a user can monitor the inputs. Also, default settings are
    set in this function and are overwritten if a user provides the parameter
    instead. 
    
    : Param program_parameters: Name of the parameter file being parsed
    
    : Return dictionary: Returns a dictionary of all paramters for later parts of the
                          program (lots of variables)
    
    """
    
    # Default parameters (will be overridden by user---input)

    # Variant Level Default Paramters
    quality_score_min= 20

    # Default parameter for indel testing
    indel_exclusion_region_length = 1

    # Sample level default parameters
    min_total_read_count=20

    # Meta-Analysis Cutoff Values
    meta_BH_adj_p_value_cutoff = 0.05
    meta_sample_p_value_cutoff = 0.05

    # Multi-Dimensional P-Value Cutoff Value
    multi_dim_adjust_pvalue_cutoff = 0.05

    # Currenting working directory default parameters (same directory as program)
    working_directory= os.getcwd()
    input_file_location = working_directory+'\\'

    # Global Output_File_Location
    output_file_location = working_directory+'/' #for UNIX environment this symbol is required and it works fine in PC submission

    # binomial probability value (50/50 Test)
    binomial_probability_value = 0.5
    
    ###LEGACY VARIABLE KEPT FOR LATER DEVELOPMENT######
    # Variables originally created to be modified, but later
    # in development realized obsolete
    min_numb_of_samples = 1

    numb_ref_alleles_allowed = 1

    numb_alt_alleles_allowed = 1


    file = open (program_parameters, 'r')

    print ("Parsed Lines from User")
    for line in file:
        if line.startswith("--"):
            line=line.rstrip('\n')
            
            print (line)
            parsed_parameters=line.split("--")
            for x in range(1, len(parsed_parameters)):
                inputs = parsed_parameters[x].split(" ")
                
                if inputs[0] == "File_Name":
                    file_name = inputs[1] 

                elif inputs[0] == "Indel_Exclusion_Region_Length":
                    indel_exclusion_region_length = inputs[1]

                    #Testing user input (verify numeric value)
                    result = test_Number_Input(indel_exclusion_region_length)
                    #If the result passes do this or if the result fails do something else
                    if result =='Pass':
                        pass
                    else:
                        print ("")
                        print ("ERROR Alert!")
                        print ("Please check Indel_Exclusion_Region_Length- incorrect input")
                        print ("Incorrect input was: ", indel_exclusion_region_length) 
                        sys.exit()

                elif inputs[0] == "Minimum_Number_of_Samples_for_ASE":
                    min_numb_of_samples = inputs[1]

                    #Testing user input (verify numeric value)
                    result = test_Number_Input(min_numb_of_samples)
                    #If the result passes do this or if the result fails do something else
                    if result =='Pass':
                        pass
                    else:
                        print ("")
                        print ("ERROR Alert!")
                        print ("Please check Minimum_Number_of_Samples_for_ASE- incorrect input")
                        print ("Incorrect input was: ", min_numb_of_samples) 
                        sys.exit()

                elif inputs[0] == "Number_of_Reference_Alleles_Allowed":
                    numb_ref_alleles_allowed = inputs[1]

                    #Testing user input (verify numeric value)
                    result = test_Number_Input(numb_ref_alleles_allowed)
                    #If the result passes do this or if the result fails do something else
                    if result =='Pass':
                        pass
                    else:
                        print ("")
                        print ("ERROR Alert!")
                        print ("Please check Number_of_Reference_Alleles_Allowed- incorrect input")
                        print ("Incorrect input was: ", numb_ref_alleles_allowed) 
                        sys.exit()

                elif inputs[0] == "Number_of_Alternative_Alleles_Allowed":
                    numb_alt_alleles_allowed = inputs[1]

                    #Testing user input (verify numeric value)
                    result = test_Number_Input(numb_alt_alleles_allowed)
                    #If the result passes do this or if the result fails do something else
                    if result =='Pass':
                        pass
                    else:
                        print ("")
                        print ("ERROR Alert!")
                        print ("Please check Number_of_Alternative_Alleles_Allowed- incorrect input")
                        print ("Incorrect input was: ", numb_alt_alleles_allowed) 
                        sys.exit()

                elif inputs[0] == "Quality_Score_Minimum_for_Variants":
                    quality_score_min = inputs[1]

                    #Testing user input (verify numeric value)
                    result = test_Number_Input(quality_score_min)
                    #If the result passes do this or if the result fails do something else
                    if result =='Pass':
                        pass
                    else:
                        print ("")
                        print ("ERROR Alert!")
                        print ("Please check Quality_Score_Minimum_for_Variants- incorrect input")
                        print ("Incorrect input was: ", quality_score_min) 
                        sys.exit()

                elif inputs[0] == "Minimum_Read_Counts":
                    min_total_read_count = inputs[1]

                    #Testing user input (verify numeric value)
                    result = test_Number_Input(min_total_read_count)
                    #If the result passes do this or if the result fails do something else
                    if result =='Pass':
                        pass
                    else:
                        print ("")
                        print ("ERROR Alert!")
                        print ("Please check Minimum_Read_Counts- incorrect input")
                        print ("Incorrect input was: ", min_total_read_count) 
                        sys.exit()

                elif inputs[0] == "Meta_BH_adj_p_value_cutoff":
                    meta_BH_adj_p_value_cutoff = inputs[1]

                    #Testing user input (verify numeric value)
                    result = test_Number_Input(meta_BH_adj_p_value_cutoff)
                    #If the result passes do this or if the result fails do something else
                    if result =='Pass':
                        pass
                    else:
                        print ("")
                        print ("ERROR Alert!")
                        print ("Please check Meta_BH_adj_p_value_cutoff- incorrect input")
                        print ("Incorrect input was: ", meta_BH_adj_p_value_cutoff) 
                        sys.exit()

                elif inputs[0] == "Meta_sample_p_value_cutoff":
                    meta_sample_p_value_cutoff = inputs[1]

                    #Testing user input (verify numeric value)
                    result = test_Number_Input(meta_sample_p_value_cutoff)
                    #If the result passes do this or if the result fails do something else
                    if result =='Pass':
                        pass
                    else:
                        print ("")
                        print ("ERROR Alert!")
                        print ("Please check Meta_sample_p_value_cutoff- incorrect input")
                        print ("Incorrect input was: ", meta_sample_p_value_cutoff) 
                        sys.exit()

                elif inputs[0] == "Multi_Dim_adjust_pvalue_cutoff":
                    multi_dim_adjust_pvalue_cutoff = inputs[1]

                    #Testing user input (verify numeric value)
                    result = test_Number_Input(multi_dim_adjust_pvalue_cutoff)
                    #If the result passes do this or if the result fails do something else
                    if result =='Pass':
                        pass
                    else:
                        print ("")
                        print ("ERROR Alert!")
                        print ("Please check FDR_Sample_p_value- incorrect input")
                        print ("Incorrect input was: ", multi_dim_adjust_pvalue_cutoff) 
                        sys.exit()

                elif inputs[0] == "Binomial_Probability_Value":
                    binomial_probability_value = inputs[1]

                    #Testing user input (verify numeric value)
                    result = test_Number_Input(binomial_probability_value)
                    #If the result passes do this or if the result fails do something else
                    if result =='Pass':
                        pass
                    else:
                        print ("")
                        print ("ERROR Alert!")
                        print ("Please check Binomial_Probability_Value- incorrect input")
                        print ("Incorrect input was: ", binomial_probability_value) 
                        sys.exit()
           
                elif inputs[0] == "Output_File_Location":
                    output_file_location = inputs[1]
    else:
        pass

    # Close the intial file
    file.close()
    
    # Printing user settings to terminal in case program crashes out before completion
    print ("")
    print ("")
    print ("Exact User Parameter Settings")
    print ("")
    print ("The input file is: ", file_name)
    print ("The output directory for analysis is: ", output_file_location)
    print ("")
    print ("")
    print ("The minimum qualtity score (phred score) for a variant is: ", quality_score_min)
    print ("The indel exclusion region length from identified indels is: ", indel_exclusion_region_length)
    #print ("The minimum number of samples to count a variant for ASE is: ", min_numb_of_samples)
    print ("The number of allowable reference alleles is (currently program is limited to one): ", numb_ref_alleles_allowed)
    print ("The number of allowable alternative alleles is (currently program is limited to one): ", numb_alt_alleles_allowed)
    print ("The minimum number of total read counts for a sample per variant is: ", min_total_read_count)
    print ("")
    print ("")
    print ("The binomial probability value for ASE testing is: ", binomial_probability_value)
    print ("")
    print ("")
    print ("Meta-Analysis of Data")
    print ("The Meta BH adjusted p-value cutoff is: ", meta_BH_adj_p_value_cutoff)
    print ("The p-value cutoff used for estimated tallying of samples is: ", meta_sample_p_value_cutoff)
    print ("")
    print ("Multi-Dimensional P-Value Adjustment")
    print ("The p-value cutoff for testing is: ", multi_dim_adjust_pvalue_cutoff)
    print ("")
    print ("")
    

    # Returns a dictionary of all the variables
    return{'file_name':file_name, 'indel_exclusion_region_length':indel_exclusion_region_length,
    'min_numb_of_samples':min_numb_of_samples, 'numb_ref_alleles_allowed':numb_ref_alleles_allowed,
    'numb_alt_alleles_allowed':numb_alt_alleles_allowed, 'quality_score_min':quality_score_min, 
    'min_total_read_count':min_total_read_count, 'binomial_probability_value': binomial_probability_value,
    'meta_BH_adj_p_value_cutoff':meta_BH_adj_p_value_cutoff, 'meta_sample_p_value_cutoff': meta_sample_p_value_cutoff,
    'multi_dim_adjust_pvalue_cutoff': multi_dim_adjust_pvalue_cutoff, 'output_file_location':output_file_location}



def parsing_input(program_parameters):

    """

    Parses apart the command line or qsub submission file to get all user input parameters
    for analyzing the data. Function also prints all user parameters to
    command line, so a user can monitor the inputs. Also, default settings are
    set in this function and are overwritten if a user provides the parameter
    instead. 

    : Param program_parameters: Name of the parameter file being parsed

    : Return dictionary: Returns a dictionary of all paramters for later parts of the
                          program (lots of variables)

    """
    
    # Default parameters (will be overridden by user---input)

    # Variant Level Default Paramters
    quality_score_min= 20

    # Default parameter for indel distance
    indel_exclusion_region_length = 1

    # Sample level default parameters
    min_total_read_count=20

    # Meta-Analysis Cutoff Values
    meta_BH_adj_p_value_cutoff = 0.05
    meta_sample_p_value_cutoff = 0.05

    # Multi-Dimensional P-Value Cutoff Value
    multi_dim_adjust_pvalue_cutoff = 0.05

    # Currenting working directory default parameters (same directory as program)
    working_directory= os.getcwd()
    input_file_location = working_directory+'\\'

    # Global Output_File_Location
    output_file_location = working_directory+'/' #for UNIX environment this symbol is required and it works fine in PC submission

    # binomial probability value (50/50 Test)
    binomial_probability_value = 0.5
    
    ###LEGACY VARIABLE KEPT FOR LATER DEVELOPMENT######
    # Variables originally created to be modified, but later
    # in development realized obsolete
    min_numb_of_samples = 1

    numb_ref_alleles_allowed = 1

    numb_alt_alleles_allowed = 1

    parsed_parameters=program_parameters.split("--")

    # print (parsed_parameters)
    for x in range(1, len(parsed_parameters)):
        inputs = parsed_parameters[x].split(" ")

        if inputs[0] == "File_Name":
            file_name = inputs[1] 

        elif inputs[0] == "Indel_Exclusion_Region_Length":
            indel_exclusion_region_length = inputs[1]

            #Testing user input (verify numeric value)
            result = test_Number_Input(indel_exclusion_region_length)
            #If the result passes do this or if the result fails do something else
            if result =='Pass':
                pass
            else:
                print ("")
                print ("ERROR Alert!")
                print ("Please check Indel_Exclusion_Region_Length- incorrect input")
                print ("Incorrect input was: ", indel_exclusion_region_length) 
                sys.exit()

        elif inputs[0] == "Minimum_Number_of_Samples_for_ASE":
            min_numb_of_samples = inputs[1]

            #Testing user input (verify numeric value)
            result = test_Number_Input(min_numb_of_samples)
            #If the result passes do this or if the result fails do something else
            if result =='Pass':
                pass
            else:
                print ("")
                print ("ERROR Alert!")
                print ("Please check Minimum_Number_of_Samples_for_ASE- incorrect input")
                print ("Incorrect input was: ", min_numb_of_samples) 
                sys.exit()

        elif inputs[0] == "Number_of_Reference_Alleles_Allowed":
            numb_ref_alleles_allowed = inputs[1]

            #Testing user input (verify numeric value)
            result = test_Number_Input(numb_ref_alleles_allowed)
            #If the result passes do this or if the result fails do something else
            if result =='Pass':
                pass
            else:
                print ("")
                print ("ERROR Alert!")
                print ("Please check Number_of_Reference_Alleles_Allowed- incorrect input")
                print ("Incorrect input was: ", numb_ref_alleles_allowed) 
                sys.exit()

        elif inputs[0] == "Number_of_Alternative_Alleles_Allowed":
            numb_alt_alleles_allowed = inputs[1]

            #Testing user input (verify numeric value)
            result = test_Number_Input(numb_alt_alleles_allowed)
            #If the result passes do this or if the result fails do something else
            if result =='Pass':
                pass
            else:
                print ("")
                print ("ERROR Alert!")
                print ("Please check Number_of_Alternative_Alleles_Allowed- incorrect input")
                print ("Incorrect input was: ", numb_alt_alleles_allowed) 
                sys.exit()

        elif inputs[0] == "Quality_Score_Minimum_for_Variants":
            quality_score_min = inputs[1]

            #Testing user input (verify numeric value)
            result = test_Number_Input(quality_score_min)
            #If the result passes do this or if the result fails do something else
            if result =='Pass':
                pass
            else:
                print ("")
                print ("ERROR Alert!")
                print ("Please check Quality_Score_Minimum_for_Variants- incorrect input")
                print ("Incorrect input was: ", quality_score_min) 
                sys.exit()

        elif inputs[0] == "Minimum_Read_Counts":
            min_total_read_count = inputs[1]

            #Testing user input (verify numeric value)
            result = test_Number_Input(min_total_read_count)
            #If the result passes do this or if the result fails do something else
            if result =='Pass':
                pass
            else:
                print ("")
                print ("ERROR Alert!")
                print ("Please check Minimum_Read_Counts- incorrect input")
                print ("Incorrect input was: ", min_total_read_count) 
                sys.exit()

        elif inputs[0] == "Meta_BH_adj_p_value_cutoff":
            meta_BH_adj_p_value_cutoff = inputs[1]

            #Testing user input (verify numeric value)
            result = test_Number_Input(meta_BH_adj_p_value_cutoff)
            #If the result passes do this or if the result fails do something else
            if result =='Pass':
                pass
            else:
                print ("")
                print ("ERROR Alert!")
                print ("Please check META_FDR_p_value- incorrect input")
                print ("Incorrect input was: ", meta_BH_adj_p_value_cutoff) 
                sys.exit()

        elif inputs[0] == "Meta_sample_p_value_cutoff":
            meta_sample_p_value_cutoff = inputs[1]

            #Testing user input (verify numeric value)
            result = test_Number_Input(meta_sample_p_value_cutoff)
            #If the result passes do this or if the result fails do something else
            if result =='Pass':
                pass
            else:
                print ("")
                print ("ERROR Alert!")
                print ("Please check META_binomial_p_value- incorrect input")
                print ("Incorrect input was: ", meta_sample_p_value_cutoff) 
                sys.exit()

        elif inputs[0] == "Multi_Dim_adjust_pvalue_cutoff":
            multi_dim_adjust_pvalue_cutoff = inputs[1]

            #Testing user input (verify numeric value)
            result = test_Number_Input(multi_dim_adjust_pvalue_cutoff)
            #If the result passes do this or if the result fails do something else
            if result =='Pass':
                pass
            else:
                print ("")
                print ("ERROR Alert!")
                print ("Please check Multi_Dim_adjust_pvalue_cutoff- incorrect input")
                print ("Incorrect input was: ", multi_dim_adjust_pvalue_cutoff) 
                sys.exit()

        elif inputs[0] == "Binomial_Probability_Value":
            binomial_probability_value = inputs[1]

            #Testing user input (verify numeric value)
            result = test_Number_Input(binomial_probability_value)
            #If the result passes do this or if the result fails do something else
            if result =='Pass':
                pass
            else:
                print ("")
                print ("ERROR Alert!")
                print ("Please check Binomial_Probability_Value- incorrect input")
                print ("Incorrect input was: ", binomial_probability_value) 
                sys.exit()
   
        elif inputs[0] == "Output_File_Location":
            output_file_location = inputs[1]

        # prints the help menu prompt
        elif inputs[0] == "help" or inputs[0]=="h" or inputs[0]=="Help":
            printing_help_menu=help_menu_prompt()

        else:
            print ("")
            print ("ERROR Alert!")
            print ("Please double check your input parameters, something was not quite right")
            print ("Type: --help to see a list of options and acceptable input for the program")
            sys.exit()
    
    # Printing user settings to terminal in case program crashes out before completion
    print ("")
    print ("")
    print ("Exact User Parameter Settings")
    print ("")
    print ("The input file is: ", file_name)
    print ("The output directory for analysis is: ", output_file_location)
    print ("")
    print ("")
    print ("The minimum qualtity score (phred score) for a variant is: ", quality_score_min)
    print ("The indel exclusion region length from identified indels is: ", indel_exclusion_region_length)
    # print ("The minimum number of samples to count a variant for ASE is: ", min_numb_of_samples)
    print ("The number of allowable reference alleles is (currently program is limited to one): ", numb_ref_alleles_allowed)
    print ("The number of allowable alternative alleles is (currently program is limited to one): ", numb_alt_alleles_allowed)
    print ("The minimum number of total read counts for a sample per variant is: ", min_total_read_count)
    print ("")
    print ("")
    print ("The binomial probability value for ASE testing is: ", binomial_probability_value)
    print ("")
    print ("")
    print ("Meta-Analysis of Data")
    print ("The Meta BH adjusted p-value cutoff is: ", meta_BH_adj_p_value_cutoff)
    print ("The p-value cutoff used for estimated tallying of samples is: ", meta_sample_p_value_cutoff)
    print ("")
    print ("Multi-Dimensional P-Value Adjustment")
    print ("The p-value cutoff for testing is: ", multi_dim_adjust_pvalue_cutoff)
    print ("")
    print ("")
    
    # Returns a dictionary of all the variables
    return{'file_name':file_name, 'indel_exclusion_region_length':indel_exclusion_region_length,
    'min_numb_of_samples':min_numb_of_samples, 'numb_ref_alleles_allowed':numb_ref_alleles_allowed,
    'numb_alt_alleles_allowed':numb_alt_alleles_allowed, 'quality_score_min':quality_score_min, 
    'min_total_read_count':min_total_read_count, 'binomial_probability_value': binomial_probability_value,
    'meta_BH_adj_p_value_cutoff':meta_BH_adj_p_value_cutoff, 'meta_sample_p_value_cutoff': meta_sample_p_value_cutoff,
    'multi_dim_adjust_pvalue_cutoff': multi_dim_adjust_pvalue_cutoff, 'output_file_location':output_file_location}
    

#########################################################################################################################
################# Universal Filtering to Identify Testable  #############################################################
#########################################################################################################################

def get_file_name(input_file_name):
    
    """

    Parses apart the pathway from the raw file name, so the file name can be passed
    to the rest of the program as a simple variable

    : Param input_file_name: Full pathway of file being analyzed
    
    : Return: The file name with pathway removed
   
    """

    input_file_name = input_file_name.split("\\")
    input_file_name = input_file_name[-1]
    input_file_name = input_file_name.split("/")
    file_name = input_file_name[-1]
    return (file_name)


def sample_information(file_name):

    """

    Function counts the number of individuals in the filtered VCF file and returns the value
    for later functions to use (Sample FDR, Meta-Analysis etc.)
    
    : Param file_name: Name of filtered file being analyzed

    : Return sample_list: List of all the samples in a file
    : Return numb_individuals: The number of individuals in the file

    """

    # Opening up variant file
    tested_variants = open(file_name, 'r')

    # Reading file line by line
    for line in tested_variants:
        
        # Getting the header information from the file
        if line.startswith('#CHROM'):
            parsed_line = line.rstrip().split('\t')

 

            # Retrieve_Sample_Names_for_Logging
            sample_list = parsed_line[9:]
            numb_individuals = len(sample_list) 

        # Getting the actual data from the file to perform meta-analysis
        else:
            break

    return {'numb_individuals': numb_individuals, 'sample_list': sample_list}


#########################################################################################################################
########################################## FDR P-value Adjustment Code ##################################################
#########################################################################################################################


"""
PROGRAM 
Implementing Benjamini-Hochberg's FDR P-value Ajustment for Multiple Testing
WRITTEN BY
M. Joseph Tomlinson
DESCRIPTION
Program implements the p-value adjustment for controlling false discovery rates
outlined in the paper by Benjamini-Hochberg (1995).

Source of code: https://github.com/mjtiv/FDR_Corrected_P-value/blob/master/FDR_Code_Matches_R
"""


def create_pvalues_fdr_results_dict(pvalues_list):


    '''

    Creates a dictionary inside a list to allow corrected
    FDR pvalues to be stored and also allow for values to changed
    and updated in the inside dictionary

    : Parameters pvalues_list: list of input pvalues

    : Return pvalues_fdr_results_list: list dictionary where the list
        value is the original p-value, followed index value then dictionary

    '''

    # Create list dictionary combo to store results
    pvalues_list_dict = []

    # Start index counter
    x = 0

    # Start looping over pvalues
    for pvalue in pvalues_list:

        # For each entry add a value as an outside list with a dictionary inside
        # Allow sorting of the outside value, very important
        pvalues_list_dict.append([pvalue, x, {'index_in_list': x,
                                               'original_pvalue': pvalue,
                                               'adjusted_pvalue': 'nan'}])
        # Add to the index counter
        x += 1

    return (pvalues_list_dict)


def reorder_sorted_pvalues_list_dict(sorted_pvalues_list_dict):

    """

    Takes in the final list-dictionary with corrected p-values
    and replaces pvalue "list" value with the original list index value
    so the list-dictionary can be re-sorted to match the original list
    order.

    : Param sorted_pvalues_list_dict: list-dictionary with all entries
        where the list value is the original p-value

    : Return re_sorted_pvalues_list_dict: list-dictionary where the
        original the list value is the original index value and has been
        re-sorted by that value

    """

    # Loop over all the entries in the list dictionary
    for entry in sorted_pvalues_list_dict:

        # Delete the first entry in the list
        # To get to the original index value
        del entry[0]

    # Sort the list (outside value using original index value)
    re_sorted_pvalues_list_dict = sorted(sorted_pvalues_list_dict)

    return (re_sorted_pvalues_list_dict)


def fdr_correction(pvalues_list):

    """

    Function performs FDR correction of the p-values
    using Benjamini-Hochberg (1995) which sorts the list of pvalues
    and then determines the p-value correction based on the rank and following
    equation (p-value x NumbTest / p-value_Rank)
    
    : Param values: list of pvalues to correct
    : Return p_value_dict: dictionary of corrected pvalues
    
    """

    # Create a list of pvalues with a dictionary inside each entry
    pvalues_list_dict = create_pvalues_fdr_results_dict(pvalues_list)

    # Sort the list in Reverse (outside value uses raw pvalue)
    # if a duplicate is found, goes to next value (aka index value),
    # index value required to prevent sort from breaking
    reverse_sorted_pvalues_list_dict = sorted(pvalues_list_dict, reverse=True)

    # Position Movement Counter
    i = 0

    # Get the total number of pvalues analyzed
    total_pvalues = len(pvalues_list_dict)

    # Start looping over list dictionary 
    for entry in reverse_sorted_pvalues_list_dict:

        # Get the sorted position of p-value (most significant to least significant),
        # opposite the current list
        p_value_sorted_position = total_pvalues - i

        # Last value in list (no need for following calculation)
        if i == 0:
            
            # Get corrected pvalue for position i
            # Important: No correction occurs for position because last value in sorted list
            fdr_adj_pvalue = round((entry[0] * total_pvalues / p_value_sorted_position), 8)

            # Adjust p-values may be greater than 1
            if fdr_adj_pvalue > 1:
                fdr_adj_pvalue = 1

            # Update the dictionary value with corrected pvalue
            entry[2]['adjusted_pvalue'] = fdr_adj_pvalue

            # Increment i value
            i+=1

        else:
            
            # Get corrected pvalue for current entry
            fdr_adj_pvalue = round((entry[0] * total_pvalues / p_value_sorted_position), 8)

            # Get pvalue for the prior entry, which has a less significant original pvalue, which after
            # adjustment could be more significant and needs to replace the adjusted pvalue
            # of the current entry if this occurs
            prior_entry_adj_pvalue = reverse_sorted_pvalues_list_dict[i -1][2]['adjusted_pvalue']

            # Correct adjusted pvalue if lower ranked pvalue is more signficant
            if prior_entry_adj_pvalue < fdr_adj_pvalue:
                fdr_adj_pvalue = prior_entry_adj_pvalue

            # Update the dictionary value with corrected pvalue
            entry[2]['adjusted_pvalue'] = fdr_adj_pvalue

            # Increment i value
            i+=1

    # Replace first value in list with index to re-order values to match orginal data
    re_sorted_pvalues_list_dict = reorder_sorted_pvalues_list_dict(reverse_sorted_pvalues_list_dict)

    return(re_sorted_pvalues_list_dict)


#########################################################################################################################
######################################## RNA-Seq Filtering Code #########################################################
#########################################################################################################################

def failure_Report(failure_file_out, parsed_line, failure):

    """

    Writes all the failure variants to a seperate file with the exact
    failure defined to allow validation of filtering
    
    : Param failure_File_Out: Failure file being written to
    : Param parsed_line: Variant data being passed to program
    : Param failure: Exact failure for the variant
    
    : Return: NONE

    """

    line=('\t'.join(map(str,parsed_line)))
    failure_file_out.write(str(failure)+"\t"+line)
    return()


def removeSNPsInExclusionZone(parsed_line, no_of_exclusions, indel_exclusion_regions):

    """

    Identifies variants that are found in INDEL exclusion zones for filtering
    
    : Param parsed_line: Variant being examined
    : Param noex: Number of INDEL exclusion zones
    : Param indel_exclusion_regions: Identified INDEL exclusion regions for data

    : Return line_of_data: variant being examined
    : Return judging_score: verdict if data is found in indel zone

    """
   
    # Judge Region Exclusion (local counter to pass or fail SNPs)
    judregex = 0
    
    # Examining the overlap between the SNP and indel regions
    for i in range(no_of_exclusions): 
   
    # Nested "if" loops to examine if SNP chromosomes match then compares the
    # the regions to see if SNP is Indel region (moves couter)
        if parsed_line[0] == indel_exclusion_regions[i][0]:
            if indel_exclusion_regions[i][1] < int(parsed_line[1]) < indel_exclusion_regions[i][2]:
                # Judge if variant is in exclusion region counter, if fail will be
                # Excluded from future analysis
                judregex += 1
                continue
    
    return{'line_of_data':parsed_line, 'judging_score':judregex}


def identify_INDEL_Regions(output_file_location, input_vcf, exclusion_region_length):
   
    """

    Identifies all INDELs in the dataset that could complicate variant analysis using information
    from both the reference and alternative alleles recorded in a vcf file
    
    : Param inputvcf: Input vcf file that is being analyzed for INDELs
    : Param exreglen: Exclusion region based on sequencing length of data
    
    : Return: A dictionary with all the exclusion regions and total number of regions excluded
    
    """

    # infile being analyzed for INDELS
    vcf_input_file = open(input_vcf, "r")

    file_name = get_file_name(input_vcf)

    indel_log_file = open(output_file_location + "/Filtering_Results/Identified_Indel_Regions.txt", "w")
    
    ###Section of Code looks at the Indels That Passed and Counts Their Numbers

    #Creating a list of exclusion regions for Indels
    indel_exclusion_regions = []
    #Start Counter for no_of_exclusions
    no_of_exclusions = 0

    # Creating a list of lengths (stats about indels)
    indel_lengths_list =[]

    for line in vcf_input_file:
        
        #Get Header from file
        if line.startswith('#CHROM'):
            parsed_line = line.split('\t')
            header_info = parsed_line[:8]
            header_info = ('\t'.join(map(str,header_info)))
            indel_log_file.write("Chromosome\tStart_INDEL_Zone\tStop_INDEL_Zone\t" + header_info + "\n")

            #indel_log_file.write("Start\tStop\t
    
        elif line.startswith(("##", "#", " #", "'#", '"##')):
                    pass
                
        else:
            parsed_line  = line.rstrip().split('\t')
            # Is looking at the column called FILTER-- if "PASS" it passed all GATK filters, if failed
            # will list the filter it failed examples DP=filtered depth issue 
            if parsed_line[6] != 'PASS':  #This filter removes TONS of samples
                continue
            else:
                #Identifying the idels in the VCF file (columns 3 and 4)
                #Retrieve the alternative allele column from VCF
                alleleR = parsed_line[3].split(',') #splitting for different types of variants
                alleleA = parsed_line[4].split(',') #splitting for different types of variants
        
                #find the INDELs in either Ref or Alt using for loops
                Indel= 0
                for i in alleleR:
                    i = i.rstrip('"')
                    i = i.lstrip('"')
                    if len(i) > 1:
                        indel_lengths_list.append(len(i))
                        Indel += 1
                #Identify alternative alleles where indels or a deletion (asterix symbol)
                for j in alleleA:
                    j = j.rstrip('"')
                    j = j.lstrip('"')
                    if len(j) > 1 or j=="*":
                        Indel += 1
                        if len(j) > 1:
                            indel_lengths_list.append(len(j))
                        else:
                            indel_lengths_list.append(-1)
                        
                #Adds the Indel to a list and calculates the foward and reverse distance of it
                if Indel > 0:
                    indel_exclusion_regions.append([])
                    indel_exclusion_regions[no_of_exclusions].append(parsed_line[0])
                    indel_exclusion_regions[no_of_exclusions].append(int(parsed_line[1]) - exclusion_region_length)
                    indel_exclusion_regions[no_of_exclusions].append(int(parsed_line[1]) + exclusion_region_length)

                    #Counter for number of exclusions
                    no_of_exclusions += 1

                    #Write Info to Log file 
                    start_exclusion_zone = int(parsed_line[1]) - exclusion_region_length
                    stop_exclusion_zone = int(parsed_line[1]) + exclusion_region_length

                    indel_log_file.write(parsed_line[0] + "\t" + str(start_exclusion_zone)
                                         + "\t" + str(stop_exclusion_zone) +"\t")

                    variant_data = (parsed_line[:8])
                    variant_data = ('\t'.join(map(str,variant_data)))

                    indel_log_file.write(variant_data + "\n")
                             
    vcf_input_file.close()  
    indel_log_file.close()

    
    # indel stats dictionary

    # Note: some variants have multiple indels, so number of indels will be different
    # from number of indel regions
    number_indels_identified = len(indel_lengths_list)

    
    if number_indels_identified > 0:
        longest_indel = max(indel_lengths_list)
        shortest_indel = min(indel_lengths_list)
        average_indel = sum(indel_lengths_list)/len(indel_lengths_list)

    # Small datasets may not have any indels (safety built in)
    else:
        longest_indel = 0
        shortest_indel = 0
        average_indel = 0


    indel_stats_dict ={}
    indel_stats_dict.update({'no_of_exclusions_regions': no_of_exclusions})
    indel_stats_dict.update({'number_indels_identified': number_indels_identified})
    indel_stats_dict.update({'longest_indel': longest_indel})
    indel_stats_dict.update({'shortest_indel': shortest_indel})
    indel_stats_dict.update({'average_indel': average_indel})
                
    return {'indel_exclusion_regions':indel_exclusion_regions, 'no_of_exclusions':no_of_exclusions,
            'indel_stats_dict': indel_stats_dict}

def testing_variant_counts(parsed_line, min_total_read_count):
    
    """

    Tests all samples for a variant to see if the variant can be futher analyzed.
    Various filtering criterias are applied to see if the variant is "testable."
    Can only examine bi-allelic samples 
    
    : Param Parsed_line: Line of data (variant) being examined
    
    : Return Dictionary: Dictionary of various results from the variant analysis
    
    """

    #Variants starts as 'fail' if it makes it through the list will become 'pass'
    verdict = 'fail'

    #Local Counters---so if all samples show behavior triggers larger counter
    homozygous_ref_variant_local_counter = 0
    homozygous_alt_variant_local_counter = 0
    low_read_count_local_counter = 0
    low_freq_count_local_counter = 0
    no_value = 0
    
    #Testable Variants Counter
    testable_SNPs = 0

    #starts the range of i to capture all samples genotype information
    #format for each sample genotype is (GT:AD:DP:GQ:PL)
    for sample in parsed_line[9:]:
        sample = sample.rstrip('"')
        sample = sample.lstrip('"')

        # Verify genotype data contains more information other files may not
        if sample.count(':'):
            
        #get the genotype information as a list for each sample (ex: ['0/1', '135,464', '605', '99', '11857,0,2154'])
            genotyping_information = sample.split(':')

            #Skipping all no recorded genotype values
            if genotyping_information[0] == './.':
                no_value += 1
                #print ("No Record")
                continue
    
            #get the alleles for a genotype (ex: ['0', '1'])
            alleles = genotyping_information[0].split('/')

            #Returning the physical genotype values
            allele_one=(alleles[0])
            allele_two=(alleles[1])

            # Dealing with empty genotyping counts values
            if genotyping_information[1] == './.' or genotyping_information[1] == '.':
                no_value += 1
                continue

            #get the count of the genotype counts (ex: ['135', '464'])
            if len(genotyping_information[1]) > 1: 
                genotyping_counts = genotyping_information[1].split(',')
            
                #Gets the allele total (ex: 135 + 464 = 599)
                total_counts = int(genotyping_counts[0])+int(genotyping_counts[1])

                #if the number of alleles is less than user input (standard value = 20) stop analysis
                if total_counts < min_total_read_count:
                    low_read_count_local_counter+= 1
                    #print ("Number of Reads Too Low")
                    continue

                #Allele count for individuals [alt, ref] (ex: 464, 135)
                count_list = []

                # Note: alleles get loaded alternate first followed by reference
                if int(genotyping_counts[0]) > 0:
                    count_list.append(int(genotyping_counts[0]))
                    
                if int(genotyping_counts[1]) > 0:
                    count_list.append(int(genotyping_counts[1]))
                    
                # Finds the minimum in the list
                lowest_allele_count = min(count_list)

                # Checks to see minimum allele is < (1% of total alleles for bi-allelic samples
                # Issue with this trigger is monoallelic samples that get one biallelic count
                if lowest_allele_count <= (0.01*total_counts) and allele_one != allele_two:
                    low_freq_count_local_counter+=1
                    #print ("Lowest Allele Count Too Low")
                    continue

            # Dealing with homozygous reads with only one reported count (occurs sometimes)
            else:
                if int(genotyping_information[1]) < min_total_read_count:    
                    low_read_count_local_counter+= 1
                    #print ("Number of Reads Too Low")
                    continue

            ####################Testing##############
            #Returning the physical genotype values
            allele_one=(alleles[0])
            allele_two=(alleles[1])

            #comparing genotypes if homozygous ref stops (ASE detection only heterzygous SNPs)
            if (allele_one == allele_two and allele_one == '0'):
                homozygous_ref_variant_local_counter+= 1
                continue

            #comparing genotypes if homozygous alt stops (ASE detection only heterzygous SNPs)
            if (allele_one == allele_two and allele_one == '1'):
                homozygous_alt_variant_local_counter+= 1
                continue
            ########################################
      
            # sample passed all thresholds---passes and can be tested
            verdict = 'pass'

        # Dealing with no data issue (like NO DATA found)
        else:
            no_value += 1
            #print ("No Record")
            continue
                                
    return {'verdict': verdict, 'homozygous_ref_variant_local_counter': homozygous_ref_variant_local_counter,
            'homozygous_alt_variant_local_counter': homozygous_alt_variant_local_counter,
            'low_read_count_local_counter': low_read_count_local_counter,
            'low_freq_count_local_counter': low_freq_count_local_counter,
            'no_value': no_value} 


def parse_filter_and_binomial_test(filtered_rna_seq_file, parsed_line, parameter_stuff):

    """

    Prints final sample data to file, but after flagging genotypes with various status
    depending on if it passes or not. All samples with "Biallelic" flag passed all the
    filters and was further analyzed for ASE.
    
    : Param filtered_rna_seq_file: Filtered file being written
    : Param paramter_stuff: Input paramter variables for testing
    : Param parsed_line: Variant data from file
    
    : Return NONE

    """

    # Input parameters for program
    min_total_read_count = int(parameter_stuff['min_total_read_count'])
    binomial_probability_value = float(parameter_stuff['binomial_probability_value'])

    variant_info = parsed_line[:8]
    variant_info = ('\t'.join(map(str,variant_info)))
    filtered_rna_seq_file.write(variant_info+"\t")
    filtered_rna_seq_file.write("Verdict:Genotype:Counts:Binomial_P_value\t")

    for sample in parsed_line[9:]:
        sample = sample.rstrip('"')
        sample = sample.lstrip('"')

        # Verify genotype data contains more information other files may not
        if sample.count(':'):

            #get the genotype information as a list for each sample (ex: ['0/1', '135,464', '605', '99', '11857,0,2154'])
            genotyping_information = sample.split(':')

            #Skipping all no recorded genotype values
            if genotyping_information[0] == './.':
                # Write to file and move on
                filtered_rna_seq_file.write("No_Data:" + str(genotyping_information[0])
                                            + ":" + str(genotyping_information[0]) + ":NA\t")
                continue

            # Dealing with empty genotyping counts values
            if genotyping_information[1] == './.' or genotyping_information[1] == '.':
                filtered_rna_seq_file.write("No_Data:" + str(genotyping_information[0])
                                            + ":" + str(genotyping_information[0]) + ":NA\t")
                continue


            #get the count of the genotype counts (ex: ['135', '464'])
            if len(genotyping_information[1]) > 1:
                
                genotyping_counts = genotyping_information[1].split(',')
            
                #get the alleles for a genotype (ex: ['0', '1'])
                alleles = genotyping_information[0].split('/')

                #Returning the physical genotype values, so can be later used to INDEX the data for
                #SNPs with 2 alternative alleles
                allele_one=(alleles[0])
                allele_two=(alleles[1])

                #get the count of the genotype counts (ex: ['135', '464'])
                genotyping_counts = genotyping_information[1].split(',')

                #Gets the allele total (ex: 135 + 464 = 599)
                total_counts = int(genotyping_counts[0])+int(genotyping_counts[1])

                # Flagging and Filtering homozygous SNPs with issues
                if allele_one == allele_two:

                    #if the number of alleles is less than user input (standard value = 20) stop analysis
                    if total_counts < min_total_read_count:
                        filtered_rna_seq_file.write("Homo_Low_Count:"+str(genotyping_information[0])
                                                    +":"+str(genotyping_information[1])+":NA\t")
                        continue

                    else:
                       filtered_rna_seq_file.write("Homo:"+str(genotyping_information[0])
                                                    +":"+str(genotyping_information[1])+":NA\t")
                       continue

                # Filtering Low Count Biallelic Variants
                if total_counts < min_total_read_count:
                        filtered_rna_seq_file.write("Low_Read_Count:"+str(genotyping_information[0])
                                                    +":"+str(genotyping_information[1])+":NA\t")
                        continue

                # Filtering and Testing Biallelic Variants
                else:
                    
                    # Allele count for individuals [alt, ref] (ex: 464, 135)
                    count_list = []

                    # Built in for safety, only really affected homozygous SNPs in other programs
                    # Note: alleles get loaded alternate first followed by reference
                    if int(genotyping_counts[0]) > 0:
                        count_list.append(int(genotyping_counts[0]))
                    if int(genotyping_counts[1]) > 0:
                        count_list.append(int(genotyping_counts[1]))
                    
                    # Finds the minimum in the list
                    lowest_allele_count = min(count_list)
                
                    # Checks to see minimum allele is < (1% of total alleles for bi-allelic samples
                    # Issue with this trigger is monoallelic samples that get one biallelic count
                    if lowest_allele_count <= (0.01*total_counts) and allele_one != allele_two:
                        filtered_rna_seq_file.write("Low_Allele_Count:" +str(genotyping_information[0])
                                                +":"+str(genotyping_information[1])+ ":NA\t")
                        continue

                    # Performs the SciPy Binominal test
                    # Reference website for test
                    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binom_test.html
                    pvalue = (stats.binom_test(int(genotyping_counts[0]),total_counts, binomial_probability_value))
                    filtered_rna_seq_file.write("Biallelic:" +str(genotyping_information[0])
                                                +":"+str(genotyping_information[1])+ ":" + str(pvalue)+ "\t")

            # Dealing with homozygous reads with only one reported count (occurs sometimes in data)
            else:
                if int(genotyping_information[1]) < min_total_read_count:    
                    filtered_rna_seq_file.write("Homo_Low_Count:"+str(genotyping_information[0])
                                                    +":"+str(genotyping_information[1])+":NA\t")
                    continue

                else:
                    filtered_rna_seq_file.write("Homo:"+str(genotyping_information[0])
                                                    +":"+str(genotyping_information[1])+":NA\t")
                    continue
                
        # Dealing with bad data that has no reported values (./.)
        else:
            sample = sample.rstrip('\n')
            filtered_rna_seq_file.write("No_Data:"+str(sample) + ":NA:NA\t")
        
    filtered_rna_seq_file.write("\n")

    return()


def filter_RNA_Seq_Data(parameter_stuff):

    """

    Filters the raw RNA-Seq data using various parameters to create a testable filtered dataset that
    furthered analyzed for ASE

    : Param input_file: Input RNA-Seq VCF file being filtered
    : Param exclusion_region_length: Length of the region surrounding indels to filter neighboring SNPs
    : Param quality_score_min: Minimum quality score to filter a variant
    : Param min_total_read_count: Minimum per sample read count

    : Return filtered_rna_seq_file_name: Name of the filtered file to use later in the program
    : Return raw_rna_seq_stats: Stats about the overall filtering process
    : Return indel_stats_dict: Statistics about the indels in the vcf file

    """

    # Get variables from parameter file
    input_file = parameter_stuff['file_name']
    indel_exclusion_region_length = int(parameter_stuff['indel_exclusion_region_length'])
    #min_numb_of_samples = int(parameter_stuff['min_numb_of_samples'])
    numb_ref_alleles_allowed = int(parameter_stuff['numb_ref_alleles_allowed'])
    numb_alt_alleles_allowed = int(parameter_stuff['numb_alt_alleles_allowed'])
    quality_score_min = int(parameter_stuff['quality_score_min'])
    min_total_read_count = int(parameter_stuff['min_total_read_count'])
    output_file_location = parameter_stuff['output_file_location']

    # Get the RNA_Seq File Name With Extra Pathway Stuff
    rna_seq_file = input_file
    # Strips pathway from file name
    rna_seq_file_name = get_file_name (rna_seq_file)

    # Open file
    rna_seq_file = open (input_file, 'r')

    #Passed Variants to be Tested
    filtered_rna_seq_file_name = output_file_location + "/Filtering_Results/Testable_Informative_Filt_Variants.txt"
    filtered_rna_seq_file = open (filtered_rna_seq_file_name, 'w')

    #Log file of variants that fail and could NOT be tested (quality controls)
    gatk_failure_file_out=open(output_file_location + "/Filtering_Results/GATK_Failing_RNA_Seq_variants.txt",'w')
    other_failure_file_out = open(output_file_location + "/Filtering_Results/Other_Failing_RNA_Seq_variants.txt",'w')

    #File Stats Dictionary
    raw_rna_seq_stats = {'variants_in_file': 0,
                        'no_of_individuals_in_file': 0,
                        'failed_GATK_SNP_filter' : 0, 'failed_Qual_filter': 0,
                        'numb_SNPs_excl_Indels' : 0,
                        'ref_allele_mutiple_forms' : 0, 'alt_allele_mutiple_forms': 0,
                        'no_genotype_values': 0, 'low_read_count': 0, 'low_freq_count': 0,
                        'combo_filter_failure': 0,
                        'no_indel_exclusion_regions': 0,
                        'no_ref_homozygous_variants' : 0,
                        'no_alt_homozygous_variants' : 0,
                        'no_combo_homozygous_variants' : 0,
                        'passing_variants': 0}
    

    #Identifies the INDEL regions in the dataset based on a certain length
    indel_data=identify_INDEL_Regions(output_file_location, input_file, indel_exclusion_region_length)
    indel_exclusion_regions = indel_data['indel_exclusion_regions']
    no_of_exclusions = indel_data['no_of_exclusions']
    indel_stats_dict = indel_data['indel_stats_dict']
    
    # Tally number of indel regions
    raw_rna_seq_stats['no_indel_exclusion_regions'] = no_of_exclusions


    #Starts parsing through the data
    for line in rna_seq_file:

        #Get the header from the file and print to new file
        if line.startswith("#CHROM"):
            file_header_info=line
            filtered_rna_seq_file.write(file_header_info)
            gatk_failure_file_out.write("Failure\t"+line)
            other_failure_file_out.write("Failure\t"+line)

            #Getting sample IDs from list
            line = line.rstrip('\n')
            line = line.split("\t")
            sample_list=line[9:]

            # Cutting the number of individuals in the file
            raw_rna_seq_stats['no_of_individuals_in_file']= len(sample_list)
            
        elif line.startswith(("##", "#", " #", "'#", '"##')):
                    pass

        else:
            #Counting the Number of Probes in File
            raw_rna_seq_stats['variants_in_file']+= 1

            #Splitting the line
            parsed_line = line.split("\t")

            #Removes all SNPs without a filter passing score
            if parsed_line[6] != 'PASS':
                raw_rna_seq_stats['failed_GATK_SNP_filter']+= 1
                failure="GATK_Filter"
                #writing GATK failures to its own file because SO MANY
                failure_Report(gatk_failure_file_out, parsed_line, failure)
                continue
            
            #Quality score filter of data (QUAL column)
            if float(parsed_line[5]) < quality_score_min:
                raw_rna_seq_stats['failed_Qual_filter']+= 1
                failure="Qual_Score_Filter"
                failure_Report(other_failure_file_out, parsed_line, failure) 
                continue

            #Filters all SNPs found in INDEL Exclusion Zone
            exclusion_filter = removeSNPsInExclusionZone(parsed_line, no_of_exclusions, indel_exclusion_regions)
            parsed_line = exclusion_filter['line_of_data']
            judregex=exclusion_filter['judging_score']
            if judregex > 0:
                raw_rna_seq_stats['numb_SNPs_excl_Indels']+= 1
                failure="Indel_Region"
                failure_Report(other_failure_file_out, parsed_line, failure)
                continue

            #Splitting the Reference Alleles
            alleles_Ref = parsed_line[3].split(',')

            #Filters if Reference Allele has more than one form by chance that INDEL filter didn't catch
            # should occur, but built in for safety
            no_Ref_Alleles = len(alleles_Ref)
            if no_Ref_Alleles > 1:
                raw_rna_seq_stats['ref_allele_mutiple_forms']+= 1
                failure="Multiple_Ref_Alleles"
                failure_Report(other_failure_file_out, parsed_line, failure)
                continue

            #Splitting the Alternative Alleles
            alleleA = parsed_line[4].split(',')

            #Filters alternative allele if greater than input threshold alleles
            noAtl = len(alleleA)
            if noAtl > 1:
                raw_rna_seq_stats['alt_allele_mutiple_forms']+= 1
                failure="To_Many_Alt_Alleles"
                failure_Report(other_failure_file_out, parsed_line, failure)
                continue
            
            # Testing if variant can actually be analyzed (one sample biallelic)
            variant_counts = testing_variant_counts(parsed_line, min_total_read_count)

            # Various filtering Variables
            homozygous_ref_variant_local_counter = variant_counts['homozygous_ref_variant_local_counter']
            homozygous_alt_variant_local_counter = variant_counts['homozygous_alt_variant_local_counter']
            low_read_count_local_counter = variant_counts['low_read_count_local_counter']
            low_freq_count_local_counter = variant_counts['low_freq_count_local_counter']
            no_value_local_counter = variant_counts['no_value']

            # True or False statement for if variant passes and should be examined further
            verdict = variant_counts['verdict']

            # Recording the type of failure of data (most SNPs fail for a variety of reasons)
            noind = len(parsed_line[9:])

            if verdict == 'fail':

                if homozygous_ref_variant_local_counter == noind:
                    raw_rna_seq_stats['no_ref_homozygous_variants'] += 1
                    failure="All Samples Ref Homozygous"
                    failure_Report(other_failure_file_out, parsed_line, failure)
                    continue

                if homozygous_alt_variant_local_counter == noind:
                    raw_rna_seq_stats['no_alt_homozygous_variants'] += 1
                    failure="All Samples Alt Homozygous"
                    failure_Report(other_failure_file_out, parsed_line, failure)
                    continue

                # Counting how many variants where failture is due to comboniation of homozygous
                # alleles
                if (homozygous_ref_variant_local_counter + homozygous_alt_variant_local_counter) == noind:
                    raw_rna_seq_stats['no_combo_homozygous_variants'] += 1
                    failure="All Samples Combo Homozygous"
                    failure_Report(other_failure_file_out, parsed_line, failure)
                    continue

                if no_value_local_counter == noind:
                    raw_rna_seq_stats['no_genotype_values'] += 1
                    failure="No_Genotype_Value"
                    failure_Report(other_failure_file_out, parsed_line, failure)
                    continue
                
                if low_read_count_local_counter == noind:
                    raw_rna_seq_stats['low_read_count'] += 1
                    failure="Samples_Low_Read_Count"
                    failure_Report(other_failure_file_out, parsed_line, failure)
                    continue

                if low_freq_count_local_counter == noind:
                    raw_rna_seq_stats['low_freq_count']+= 1
                    failure="Samples_Low_Freq"
                    failure_Report(other_failure_file_out, parsed_line, failure)
                    continue
            
                if (homozygous_ref_variant_local_counter
                    + homozygous_alt_variant_local_counter
                    + low_read_count_local_counter
                    + low_freq_count_local_counter + no_value_local_counter)== noind:
                    raw_rna_seq_stats['combo_filter_failure']+= 1    
                    failure="Sample_Combo_Filter"
                    failure_Report(other_failure_file_out, parsed_line, failure)
                    continue

            # Final Results Print to File for further testing
            else:
                raw_rna_seq_stats['passing_variants']+= 1
                parse_filter_and_binomial_test(filtered_rna_seq_file, parsed_line, parameter_stuff)
  

    #Flusing Data if writing slow down
    rna_seq_file.flush()
    
    #Closing Files            
    rna_seq_file.close()
    filtered_rna_seq_file.close()
    gatk_failure_file_out.close()
    other_failure_file_out.close()

    return {'filtered_rna_seq_file_name': filtered_rna_seq_file_name, 'raw_rna_seq_stats':raw_rna_seq_stats, 'indel_stats_dict': indel_stats_dict}
    

#########################################################################################################################
#################################### Sample Multi-dimensional Pvalue Adjustment #########################################
#########################################################################################################################

def determine_passing_pvalues(lowest_pvalue_list,
                              fdr_corrected_dict,
                              multi_dim_adjust_pvalue_cutoff):

    """

    Loops over the dictionary of pvalues using the original
    uncorrected pvalue list to control order and determines
    all FDR corrected pvalues that pass the cutoff

    : Param lowest_pvalue_list: list of lowest pvalues
    : Param fdr_corrected_dict: dictionary of all corrected pvalues
    : Param multi_dim_adjust_pvalue_cutoff: cutoff for determining p-value significance

    : Return count_passing_pvalues: count of all pvalues that pass the cutoff

    """

    # List to store all passing pvalues
    passing_pvalues = []

    # Inex in List Counter for Movement
    index_in_list_counter = 0

    # Looping over the dictionary
    for value in lowest_pvalue_list:

        # Get the dictionary results part of the list-dictionary combo
        dictionary_results =(fdr_corrected_dict[index_in_list_counter][1])

        # Add to index_in_list_counter now that values have been retrieved
        index_in_list_counter += 1

        # Get the original pvalue for double-checking it matches
        original_pvalue = dictionary_results['original_pvalue']

        # Get the adjusted pvalue
        fdr_pvalue = dictionary_results['adjusted_pvalue']

        # Safety to kill program if lists being sorted over index order is wrong
        if original_pvalue != value:
            print ("P-value Lists Order is not correct, double check code")
            print ("Killing program")
            sys.exit(1)

        # Determine if pvalue is below the 0.05 cutoff
        if fdr_pvalue < float(multi_dim_adjust_pvalue_cutoff):

            # Add to the list if it passes
            passing_pvalues.append(fdr_pvalue)

        else:
            continue

    # Get the count of passing values
    count_passing_pvalues = len(passing_pvalues)

    return (count_passing_pvalues)


def determine_passing_FDR_pvalues(input_file_name, multi_dim_adjust_pvalue_cutoff):

    """

    Opens the inputted file and reads through the data extracting
    the lowest pvalue for each variant and Bonferroni corrects. The final
    pvalue is then FDR corrected and then the final list of FDR corrected
    pvalues examined to see how many pass 0.05 cutoff.

    : Param input_file_name: name of input file to open and analyze
    : Param multi_dim_adjust_pvalue_cutoff: cutoff for determining p-value significance

    : Return passing_count: number of variants that pass overall correction
    : Return variant_counter: number of variants analyzed

    """

    # Open the file
    input_file = open(input_file_name, 'r')

    # Line movement counter
    line_counter = 0

    # Lowest p-value list
    lowest_bonf_pvalue_list = []

    # Total variants count
    variant_counter = 0

    # Read through the lines of the file
    for line in input_file:

        # Skip over the header line
        if line_counter == 0:
            line_counter += 1
            continue

        # Start parsing variant results
        else:

            # Start counting variants
            variant_counter +=1

            # Clean up line
            line = line.rstrip("\n")

            # Split the data
            data = line.split("\t")

            # Get all the samples data 
            samples_data = data[9:]

            # Variant Information 
            variant_info = data[:9]

            # Create an empty list to store pvalues
            samples_pvalues = []

            # Start looping over actual data
            for sample in samples_data:

                # Clean up sample data some (remove quote marks)
                sample = sample.lstrip('"')
                sample = sample.rstrip('"')

                # Filter out all samples that are not Biallelic
                if not sample.startswith("Biallelic"):
                    continue

                # Split the sample data up to get pvalues
                sample_data = sample.split(":")

                # Get the pvalue from the results list
                sample_pvalue = float(sample_data[3])

                # Add pvalue to sample pvalues list for variant
                samples_pvalues.append(sample_pvalue)

            # Find lowest p-value
            min_pvalue = np.amin(samples_pvalues)

            # Length of list
            length_of_list = len(samples_pvalues)

            # Calculate the bonferroni pooled p-value
            bonf_pooled_pvalue = min_pvalue * length_of_list

            # Adjust bonferroni value if > 1
            if bonf_pooled_pvalue > 1:
                bonf_pooled_pvalue = 1
                
            # Add value to list for FDR correction
            lowest_bonf_pvalue_list.append(bonf_pooled_pvalue)

    # Close the input file
    input_file.close()

    # Get FDR corrected p-values
    fdr_corrected_dict = fdr_correction(lowest_bonf_pvalue_list)

    # Determine the number of passing p-values after FDR correction
    # Passes the multi_dim_adjst_pvalue_cutoff to function here
    passing_count = determine_passing_pvalues(lowest_bonf_pvalue_list,
                                              fdr_corrected_dict,
                                              multi_dim_adjust_pvalue_cutoff)
                                      
    return(passing_count, variant_counter)


def determine_biallelic_samples(samples_data):

    """

    Determine the number of biallelic samples for a variant

    : Parameter samples_data: sample data for line

    : Return biallelic_sample_line_count: number of biallelic samples

    """

    # Counter for biallelic samples
    biallelic_sample_line_count = 0

    # Start looping over values to identify biallelic samples
    for sample in samples_data:

        # Clean up sample data some (remove quote marks)
        sample = sample.lstrip('"')
        sample = sample.rstrip('"')

        # Identify biallelic samples
        if sample.startswith("Biallelic"):

            # Add to the counter
            biallelic_sample_line_count += 1

    return (biallelic_sample_line_count)
    
        
def analyze_variants_for_significance_multi_dim_test(input_file_name, passing_variant_count,
                                                     total_variants_analyzed, multi_dim_sample_output_dir,
                                                     multi_dim_adjust_pvalue_cutoff):

    """

    Parses through the input file and calculates the multi-dimensional significance threshold
    for each variant and then examines each biallelic sample for significance for that
    variant and returns pass (sig) or fail (not_sig) for each p-value in a final report with all prior
    information in the file.

    : Param input_file_name: name of the file being analyzed
    : Param passing_variant_count: number of variants that pass p-value cutoff for significance
    : Param total_variants_analyzed: total number of variants that are found in the file
    : Param multi_dim_sample_output_dir: directory where to write the results
    : Param multi_dim_adjust_pvalue_cutoff: cutoff for determining p-value significance

    : Return output_file_name: name of output file from analysis


    """

    # Open the input file
    input_file = open(input_file_name, 'r')

    # Output file name
    output_file_name = multi_dim_sample_output_dir + "/multi_dim_adj_pvalues_testable.txt"

    # Create an output file for results
    output_file = open(output_file_name, 'w')

    # Start looping over data
    for line in input_file:

        # Identify the header line
        if line.startswith("#CHROM"):

            # Write header to the new file
            output_file.write(line)

        # Loop over rest of data
        else:

            # Clean up line
            line = line.rstrip("\n")
            line = line.rstrip("\t")

            # Split the data
            data = line.split("\t")

            # Variant Information 
            variant_info = data[:8]

            # Get format information
            format_info = data[8]

            # Get all the samples data 
            samples_data = data[9:]

            # Determine the number of biallelic testable samples for variant
            biallelic_sample_line_count = determine_biallelic_samples(samples_data)

            # Calculate Significance Threshold for variant
            variant_sig_threshold = round(((passing_variant_count * float(multi_dim_adjust_pvalue_cutoff)) / (biallelic_sample_line_count * total_variants_analyzed)),8)

            # Convert variant info to string
            variant_info_str = "\t".join(map(str, variant_info))

            # Print the variant info to file
            output_file.write(variant_info_str + "\t")

            # Clean up the format info
            format_info = format_info.rstrip("\t")

            # Add to the format a new input (Variant Significance Threshold)
            new_format_info = (format_info + ":Var_Sig_Thres<" + str(variant_sig_threshold) + "\t")

            # Write new format information to file
            output_file.write(new_format_info)

            # Start looping over the sample data using new threshold
            for sample in samples_data:

                # Clean up sample data some (remove quote marks)
                sample = sample.lstrip('"')
                sample = sample.rstrip('"')

                # Identify Biallelic Samples
                if sample.startswith("Biallelic"):

                    # Split the results up
                    split_sample_results = sample.split(":")

                    # Examine if Sample less than new significance threshold
                    if float(split_sample_results[3]) < variant_sig_threshold:

                        # Remove tab before adding
                        sample = sample.rstrip("\t")

                        # Output new information
                        output_file.write(sample + ":Pass\t")

                    # If results are not significant
                    else:
                        # Remove tab before adding
                        sample = sample.rstrip("\t")

                        # Output new information
                        output_file.write(sample + ":Fail\t")
                        
                # Add NA to non-biallelic samples
                else:
                    # Remove tab before adding
                    sample = sample.rstrip("\t")

                    # Output new information
                    output_file.write(sample + ":Na\t")

            # Write new line to the file
            output_file.write("\n")

    # Nothing to return at this time
    return(output_file_name)


def filter_for_multi_dim_sig_samples(multi_dim_sample_output_dir, multi_dim_adj_file_name):

    """

    Function goes through the multi-dimensional adjusted samples and pulls out all variants with
    at least one sample showing significant ASE based on the multi-dimensional p-value adjustment
    and prints all the results to two files: variants with significant samples and variants with
    non-significant samples

    : Param multi_dim_adj_file_name: Name of the multi dimensional corrected file
    : Param multi_dim_sample_output_dir: location where to output the results

    : Return sig_multi_dimensional_file_name: file of all significant results based on adjustment
    : Return multi_dim_global_counts_dict: Tallies for global reporting

    """

    multi_dim_global_counts_dict = {
        'Biallelic_Testable_Variants': 0,
        'Total_Tests': 0,
        'Total_Sig_ASE_Variants': 0,
        'Total_Biallelic_Samples' : 0,
        'Total_Sig_Biallelic_Samples' : 0,
        'Total_Sig_ASE_Ref': 0,
        'Total_Sig_ASE_Alt': 0,
        'Total_Biallelic_No_ASE': 0,
        'Total_Passing_Homozygous': 0,
        'Total_Passing_Homozygous_Ref': 0,
        'Total_Passing_Homozygous_Alt': 0,
        'Total_Non_Testable': 0} 

    # Create a new file of only significant samples
    sig_multi_dimensional_file_name = multi_dim_sample_output_dir + "/sig_multi_dim_adj_results.txt"

    

    # Create output folder for failing data for verification purposes    
    not_sig_multi_dimensional_file_name = multi_dim_sample_output_dir + "/not_sig_multi_dim_adj_results.txt"

    # Open the files
    sig_multi_dimensional_file = open(sig_multi_dimensional_file_name, 'w')
    not_sig_multi_dimensional_file = open(not_sig_multi_dimensional_file_name,'w')
    
    # Open up the annotated multi-dimensional p-value adjusted file
    multi_dim_file = open(multi_dim_adj_file_name, 'r')
    
    #reading file line by line
    for line in multi_dim_file:

        # Setup a trigger for reporting only variants with significant sampels
        verdict = 'noASE'

        #Getting the header information from the file
        if line.startswith('#CHROM'):
    
            #Print header of file to new data_file
            sig_multi_dimensional_file.write(line)
            not_sig_multi_dimensional_file.write(line)
            continue
              
        else:
            #Splitting the data
            parsed_line = line.rstrip().split('\t')

            #Getting Variant Results Only
            variant_sample_results=parsed_line[9:]

            #Starting to Parse Sample Data
            # Counting variants in file, all should have one biallelic sample
            # from prior filtering
            multi_dim_global_counts_dict['Biallelic_Testable_Variants'] += 1
            
            for sample_data in variant_sample_results:
                split_sample_data = sample_data.split(':')
                #print (split_sample_data)

                if split_sample_data[0]=='Biallelic':
                    sample_sig_verdict = split_sample_data[4]
                    if sample_sig_verdict == 'Pass':
                        verdict = 'ASE_hit'
                    else:
                        verdic = 'No_ASE_hit'
                        pass
                    
                #not biallelic pass for now
                else:
                    pass

        #Print line of file if one ASE sample found
        if verdict == 'ASE_hit':
            multi_dim_global_counts_dict ['Total_Sig_ASE_Variants'] += 1
            sig_multi_dimensional_file.write(line)
            

        if verdict == 'noASE':
            not_sig_multi_dimensional_file.write(line)
            
    # Closing files
    multi_dim_file.close()
    sig_multi_dimensional_file.close()
    not_sig_multi_dimensional_file.close()

    return {'sig_multi_dimensional_file_name': sig_multi_dimensional_file_name,
            'multi_dim_global_counts_dict': multi_dim_global_counts_dict}


def get_list_sig_variants(file_name):

    """

    Parses through the significant variants file
    and creates a dictionary of variants and their
    information for later parsing
    
    : Param file_name: Name of file being parsed

    : Return sig_variants_dict: A dictionary of significant results  

    """

    sig_variants_dict = {}

    sig_variants_file = open(file_name, 'r')

    # Skipping the header
    for line in sig_variants_file:
        if line.startswith('#CHROM'):
            continue
        
        # Getting the actual data to create a dictionary
        else:
            parsed_line = line.rstrip().split('\t')

            variant_name = parsed_line[0] + ":" + parsed_line[1]
            variant_info = parsed_line[:8]
            sig_variants_dict.update({variant_name: variant_info})

    # Closing the file
    sig_variants_file.close()

    return(sig_variants_dict)


def multi_dim_data_for_mapping_bias_plots(testable_variants_file, sig_variants_dict, output_directory):

    """
    Function parses through all the testable biallelic variants
    and tallies the overall read counts to investigate reference allele bias

    : Param testable_variants_file: File with ALL biallelic testable variants
    : Param output_directory: Tells where to output the data file

    : Return avg_ref_allele_ratio: Average reference allele ratio ---help identify ref allele bias in data

    """

    #opening up variant file
    testable_variants_file=open(testable_variants_file, 'r') 

    #Creating a log file for variant meta analysis
    ref_bias_file=open(output_directory + "/data_for_ref_allele_bias_plotting.txt",'w')

    # Writing the heading to file for plotting file
    ref_bias_file.write("Variant_Name\tStatus\tRef_Count\tAlt_Count\tTotal_Count\tRatio\n")

    # Creating Ref Allele Ratio List
    ref_allele_ratio_list = []
    
    #reading file line by line
    for line in testable_variants_file:

        #Getting the header information from the file
        if line.startswith('#CHROM'):
            continue

        else:
            #Splitting the data
            parsed_line = line.rstrip().split('\t')

            variant_name = parsed_line[0] + ":" + parsed_line[1]

            # Determining if the variant showed ASE
            if variant_name in sig_variants_dict:
                status = 'Sig_ASE'
            else:
                status = 'No_ASE'

            # Simple Counters
            variant_total_reference_count = 0
            variant_total_alternative_count = 0
            variant_total_count = 0 

            #Getting Variant Results Only
            variant_sample_results=parsed_line[9:]

            #Starting reading and writing 
            for sample_data in variant_sample_results:
                split_sample_data = sample_data.split(':')

                if split_sample_data[0]=='Biallelic':

                    counts = split_sample_data[2].split(',')

                    reference_count = int(counts[0])
                    alternative_count = int(counts[1])
                    total_count = reference_count + alternative_count

                    #add to variant totals
                    variant_total_reference_count += reference_count
                    variant_total_alternative_count += alternative_count
                    variant_total_count += total_count 
 
                else:
                    continue

            # Get the reference allele ratio
            ref_allele_ratio = (variant_total_reference_count/variant_total_count)

            ref_allele_ratio_list.append(ref_allele_ratio)

            # Printing alll results to file
            ref_bias_file.write(str(variant_name) + "\t" + status + "\t" + str(variant_total_reference_count)
                                + "\t" + str(variant_total_alternative_count)
                                + "\t" + str(variant_total_count) + "\t"
                                + str(ref_allele_ratio) + "\n" )

    #Getting the average ref allele ratio
    avg_ref_allele_ratio = np.mean(ref_allele_ratio_list)

    # Closing the file
    testable_variants_file.close()
    ref_bias_file.close()

    return(avg_ref_allele_ratio)


def tally_final_multi_dim_results(sig_multi_dimensional_file_name, multi_dim_global_counts_dict):

    """

    Function parses through and analyzes all the results from the Multi-dimensional p-value adjustment
    and tallies the final results for final printing of reports

    : Param file_name: Name of file being analyzed (made generic so printing function can be used by other analysis)
    : Param tally_global_dictionary: Global tallying dictionary from the analysis process
    
    : Return samples_list: List of samples
    : Return samples_counters_dict: Dictionary with information about sample results
    : Return variant_list: List of variants being analyzed
    : Return variant_results_dictionary: Tallying dictionary for variants
    : Return variant_header: Header of file with variants
    : Return tally_global_dictionary: Global tallying dictionary filled in with values

    """

    # Rename Global dictionary
    tally_global_dictionary = multi_dim_global_counts_dict

    #opening up variant file
    final_results_file = open(sig_multi_dimensional_file_name, 'r')

    #Create a variant dictionary
    variant_results_dictionary = {}

    #Create a variant list for cycling through
    variant_list = []
    
    #reading file line by line
    for line in final_results_file:

        #Getting the header information from the file
        if line.startswith('#CHROM'):
            
            #Splitting the data
            parsed_line = line.rstrip().split('\t')

            # Getting Variant Header
            variant_header = parsed_line[:8]
            variant_header = ('\t'.join(map(str,variant_header)))

            # Retrieve_Sample_Names_for_Logging
            samples_list = parsed_line[9:]
            numb_individuals = len(samples_list)

            #Create a tallying dictionary for results
            samples_counters_dict = create_sample_tallying_counters(samples_list)

            pass
   
        #Getting Variant Results Only
        else:
            
            #Splitting the data
            parsed_line = line.rstrip().split('\t')
            
            #Variant information
            variant_info = parsed_line[:8]
            variant_key = ('\t'.join(map(str,variant_info)))

            # Add new variant to list
            variant_list.append(variant_key)

            variant_results_dictionary.update({variant_key: {'Biallelic_Testable': 0, 'Sig_ASE': 0, 
                           'Sig_ASE_Ref': 0, 'Sig_ASE_Alt': 0,
                           'Biallelic_No_ASE': 0,
                           'Passing_Homozygous': 0,
                           'Homozygous_Ref': 0,
                           'Homozygous_Alt': 0,                                  
                           'Non_Testable': 0}}) 
            
            # Getting variant results
            variant_sample_results=parsed_line[9:]

            #Starting reading and writing 
            x=0
            for sample_data in variant_sample_results:

                #get sample name from Index Location
                sample = samples_list[x]

                # Tallying all tests possible
                tally_global_dictionary ['Total_Tests'] += 1

                #Splitting the data for later parsing
                split_sample_data = sample_data.split(':')

                if split_sample_data[0]=='Biallelic':
                    tally_global_dictionary ['Total_Biallelic_Samples'] += 1
                    samples_counters_dict[sample]['Biallelic_Testable'] += 1
                    variant_results_dictionary[variant_key]['Biallelic_Testable'] += 1

                    # Getting Sample Verdict for Multi-Dimensional Adjustment
                    sample_sig_verdict = split_sample_data[4]

                    #Test if breach p_value cutoff
                    if sample_sig_verdict == 'Pass':
                        tally_global_dictionary ['Total_Sig_Biallelic_Samples'] += 1
                        samples_counters_dict[sample]['Sig_ASE'] += 1
                        variant_results_dictionary[variant_key]['Sig_ASE'] += 1

                        #Look to see if ref or alt is higher
                        variant_counts = (split_sample_data[2])
                        split_counts = variant_counts.split(',')

                        # If Ref is Higher Than Alt
                        if int(split_counts[0]) > int(split_counts[1]):
                            samples_counters_dict[sample]['Sig_ASE_Ref'] += 1
                            tally_global_dictionary ['Total_Sig_ASE_Ref'] += 1
                            variant_results_dictionary[variant_key]['Sig_ASE_Ref'] += 1

                            # Move Counter
                            x+= 1

                        # If Alt is Higher Than Ref
                        else:
                            samples_counters_dict[sample]['Sig_ASE_Alt'] += 1
                            tally_global_dictionary ['Total_Sig_ASE_Alt'] += 1
                            variant_results_dictionary[variant_key]['Sig_ASE_Alt'] += 1

                            # Move Counter
                            x+= 1
                            
                    #Biallic but not significant (fails cutoff)
                    else:
                        samples_counters_dict[sample]['Biallelic_No_ASE'] += 1
                        tally_global_dictionary ['Total_Biallelic_No_ASE'] += 1
                        variant_results_dictionary[variant_key]['Biallelic_No_ASE'] += 1

                        # Move Counter
                        x+= 1
                          
                #Homozygous passing sample 
                elif split_sample_data[0]=='Homo':
                    samples_counters_dict[sample]['Passing_Homozygous'] += 1
                    tally_global_dictionary ['Total_Passing_Homozygous'] += 1
                    variant_results_dictionary[variant_key]['Passing_Homozygous'] += 1

                    #Tally Type of Homozygous Genotype (allow freq calculations) 
                    homozygous_genotype = (split_sample_data[1])

                    if homozygous_genotype == "0/0":
                        samples_counters_dict[sample]['Passing_Homozygous_Ref'] += 1
                        variant_results_dictionary[variant_key]['Homozygous_Ref'] += 1
                        tally_global_dictionary ['Total_Passing_Homozygous_Ref'] += 1
                        # Move Counter
                        x+= 1

                    elif homozygous_genotype == "1/1":
                        samples_counters_dict[sample]['Passing_Homozygous_Alt'] += 1
                        variant_results_dictionary[variant_key]['Homozygous_Alt'] += 1
                        tally_global_dictionary ['Total_Passing_Homozygous_Alt'] += 1
                        # Move Counter
                        x+= 1

                    else:
                        print ("Danger Incorrect Genotype")
                        print ("Printing line of problematic data")
                        print (line)
                        print ("Printing problematic genotype")
                        print (homozygous_genotype)
                        sys.exit()
                            
                # No passing sample  
                else:
                    samples_counters_dict[sample]['Non_Testable'] += 1
                    tally_global_dictionary['Total_Non_Testable'] += 1
                    variant_results_dictionary[variant_key]['Non_Testable'] += 1

                    # Move Counter
                    x+= 1

    # Closing the initial file
    final_results_file.close()

    return (samples_list, samples_counters_dict, variant_list, variant_results_dictionary, variant_header,
            tally_global_dictionary)



#########################################################################################################################
###################################### Meta Analysis of Data  ###########################################################
#########################################################################################################################


def meta_analysis(meta_analysis_output_folder, testable_variants_file):
    
    """

    Function performs Fisher's Method to combine p-values from ASE testing
    the function first parses through the ASE log file to find all SNPs classified
    as ASE and then combines their p-values and testing using the
    scip.stats.combine_pvalues module
    
    : Param meta_analysis_output_folder: Location to put output files
    : Param testable_variants_file: Variant file with all biallelic testable variants
    
    : Return meta_analysis_all_p_values_results: Meta-analysis p-values for FDR correction
    : Return intermediate_meta_log_file_name: Name of the intermediate log file that will be re-written
        with FDR corrected q-values and then deleted
    """

    # opening up variant file
    testable_variants_file=open(testable_variants_file, 'r')

    # meta log file name
    intermediate_meta_log_file_name = meta_analysis_output_folder + '/intermediate_variants_meta_analysis_log_file.txt'

    #Creating an intermediate log file for variant meta analysis (will delete)
    inter_meta_log_file = open(intermediate_meta_log_file_name,'w')

    # create an empty list to put all result values used in analysis
    meta_analysis_all_p_values_results = []

    # reading file line by line
    for line in testable_variants_file:

        # create an empty list to put all values used in analysis
        pvalues_being_tested=[]

        # Log of samples and values tested, printed into log file
        samples_values_tested=[]
        
        
        # Getting the header information from the file
        if line.startswith('#CHROM'):

            # parsing the header line of the file
            parsed_line = line.rstrip().split('\t')

            # Get the important variant info (discard column "Mod_Format")
            file_header = parsed_line[0:7]
            
            file_header = ('\t'.join(map(str,file_header)))

            inter_meta_log_file.write(file_header+
                           '\tTotal_Samples\tAnalyzed_Values\tBi-Allelic_Samples\tchi-square_value'+
                           '\tDegrees_of_Freedom(2n)\tMeta_p-value\n')

            # Retrieve_Sample_Names_for_Logging
            sample_names = parsed_line[9:]

        #Getting the actual data from the file to perform meta-analysis
        else:
            
            # start parsing the line
            parsed_line=line.rstrip().split('\t')

            # Getting variant info
            variant_info=parsed_line[0:7]
            variant_info=('\t'.join(map(str,variant_info)))

            # Writing variant info to log file
            inter_meta_log_file.write((variant_info+('\t')))

            # Parsing apart data for sample results
            variant_results = parsed_line[9:]
            
            for x in range(len(variant_results)):

                    # Finding all testable variants (aka hetereozygous)
                    if variant_results[x].startswith("Biallelic"):
                        data=variant_results[x].split(":") 
                        p_value=data[3]

                        # getting the sample analyze from sample list
                        sample = sample_names[x]
                        
                        # Converts the string p-value to a float for fisher_testing
                        pvalues_being_tested.append(float(p_value))

                        # Create a log list for validation, round numeric values for simplicty
                        log_sample_p_value=sample+":"+str(round(float(p_value),5))

                        # Adding values examined to log
                        samples_values_tested.append(log_sample_p_value)

                    # just skipping all other types of variants in samples that fail QC
                    else:
                        continue

            # Record number of samples overall
            inter_meta_log_file.write(str(len(variant_results))+'\t')
                    
            # Creating log output for file (what samples meta-tested)
            values_tested=(', '.join(map(str,samples_values_tested)))
            inter_meta_log_file.write(values_tested+'\t')

            # Counting the number of samples tested in meta-analysis
            numb_samples=len(pvalues_being_tested)
            inter_meta_log_file.write(str(numb_samples)+'\t')

            #Performing analysis of data
            #Results are returned as (chi-square value, p-value)
            # Link for Scipy website: https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.combine_pvalues.html
            meta_analysis_results = combine_pvalues(pvalues_being_tested)

            #Printing the chi-square results
            inter_meta_log_file.write(str(meta_analysis_results[0])+'\t')

            #Printing the degrees of freedom to the file for chi-square test
            #Degrees of freedom calculated (Number of Samples X 2)
            inter_meta_log_file.write(str(numb_samples*2)+'\t')

            
            # Get the variant result p-value from meta-analysis
            meta_analysis_variant_result = (meta_analysis_results[1])

            # Convert the meta-results to 12 floating points (match normal print output from python)
            # Conversion needs to be done, otherwise dictionary and output files do not match in p-values
            # for dictionary look-up
            meta_analysis_variant_result = format(meta_analysis_variant_result, '.12f')

            # Print p-value to file
            inter_meta_log_file.write(str(meta_analysis_variant_result)+'\n')

            # Record the list of all p-values for FDR correction
            meta_analysis_all_p_values_results.append(float(meta_analysis_variant_result))

    # Flushing the data for safety reasons
    inter_meta_log_file.flush()

    # Closing meta_Log
    inter_meta_log_file.close()

    # Return values for tallying
    return {'meta_analysis_all_p_values_results': meta_analysis_all_p_values_results,
            'intermediate_meta_log_file_name': intermediate_meta_log_file_name}


def meta_printing_q_values(meta_analysis_output_folder, meta_analysis_results, p_value_threshold):

    """

    Function goes through data and prints out the corrected q-value from the FDR
    analysis

    : Param meta_analysis_output_folder: Location of the meta-analsysis results
    : Param meta_analysis_results: Results from the meta-analysis (dictionary of multiple variables)
    : Param p_value_threshold: Meta-analysis FDR p-value for significance or not

    : Return meta_analysis_sig_variants: List of significant variants from analysis
    : Return variant_meta_analysis_results_file_name: File name of final meta-anlysis folder with all output

    """
    meta_variants_results_dict = {}

    # Get the results from the meta analysis of the data
    meta_analysis_all_p_values_results = meta_analysis_results['meta_analysis_all_p_values_results']
    intermediate_meta_log_file_name = meta_analysis_results['intermediate_meta_log_file_name']

    # FDR correct the meta data
    corrected_p_values_dict = fdr_correction(meta_analysis_all_p_values_results)

    # Variants with significant values list (used for later filtering)
    meta_analysis_sig_variants_list=[]

    # open up the temporary file for analysis (reads data)
    initial_meta_log_file = open(intermediate_meta_log_file_name,'r')

    variant_meta_analysis_results_file_name = meta_analysis_output_folder + '/variant_meta_analysis_results.txt'

    # Final file for writing result to
    final_meta_analysis_file = open(variant_meta_analysis_results_file_name,'w')

    # Index in List Counter for Movement
    index_in_list_counter = 0

    # Parsing through the data
    for line in initial_meta_log_file:

        #Getting the header information from the file
        if line.startswith('#CHROM'):

            #Get the line and remove newline character
            header_line=line.rstrip() 

            final_meta_analysis_file.write(header_line+"\tBH_Adj_pvalue\tSignificant\n")

        #go through the rest of the lines
        else:

            #Removing the new line character
            line = line.rstrip('\n')

            #parse the data line
            parsed_line=line.split('\t')

            #variant name (used for tallying significant variants)
            variant=(parsed_line[0])+":"+(parsed_line[1])

            # Retrieving variants analyzed in meta-analysis for logging
            analyzed_values = parsed_line[8]

            # Retrieving chi-square value from meta-analysis (important inf values)
            chi_square_value = parsed_line[10]

            # Retrieving degrees of freedom from meta-analysis
            degrees_of_freedom = parsed_line[11]
            
            #Retrieving the original value to look up in dictionary
            p_value = (parsed_line[12]).rstrip('\t')
            p_value = float(p_value)
 
            # Get the dictionary results part of the list-dictionary combo
            dictionary_results =(corrected_p_values_dict[index_in_list_counter][1])

            # Add to index_in_list_counter now that values have been retrieved
            index_in_list_counter += 1
        
            # Get the original pvalue for double-checking it matches
            original_pvalue = dictionary_results['original_pvalue']

            # Get the adjusted pvalue
            fdr_pvalue = dictionary_results['adjusted_pvalue']

            # Safety to kill program if lists being sorted over index order is wrong
            if original_pvalue != p_value:
                print ("P-value Lists Order is not correct, double check code")
                print ("Killing program")
                sys.exit(1)

            #Rename Variable for Rest of Code
            q_value = fdr_pvalue

            final_meta_analysis_file.write(line + "\t" + str(q_value) + "\t")

            #Examine Significance for counting
            if float(q_value) < float(p_value_threshold):
                verdict = 'yes'
                final_meta_analysis_file.write(verdict + "\n")

                #get list of sig variants
                meta_analysis_sig_variants_list.append(variant)

            else:
                verdict = 'no'
                final_meta_analysis_file.write(verdict + "\n")

            meta_variants_results_dict.update({variant: {
                'chi_square_value': chi_square_value,
                'degrees_of_freedom': degrees_of_freedom,
                'analyzed_values': analyzed_values,
                'original_p_value': p_value,
                'q_value': q_value,
                'verdict': verdict}})

    # Close the intermediate file
    initial_meta_log_file.close()

    # Delete the intermediate file
    os.remove(intermediate_meta_log_file_name)

    #Close the final meta_anlaysis file
    final_meta_analysis_file.close()

    #Return the significant SNPs Counter
    return{'meta_analysis_sig_variants_list': meta_analysis_sig_variants_list,
           'variant_meta_analysis_results_file_name': variant_meta_analysis_results_file_name,
           'meta_variants_results_dict': meta_variants_results_dict}


def filter_for_meta_sig_variants(testable_variants_file_name, meta_analysis_output_folder,
                                 final_meta_analysis_sig_variants_results):

    """

    Function goes through the testable  samples and pulls out all variants that were
    found to be significant from the meta analysis and parses the testable variants
    file into two types, significant variants and non-significant variants

    : Param testable_variants_file_name: name of file with all testable variants
    : Param meta_analysis_output_folder: location where to output the data
    : Param final_meta_analysis_sig_variants_results: list of all significant variants

    : Return sig_meta_file_name: file of all significant variants (original data)
    : Return meta_global_counts_dict: Tallies for global reporting

    """

    meta_global_counts_dict = {
        'Biallelic_Testable_Variants': 0,
        'Total_Tests': 0,
        'Total_Sig_ASE_Variants': 0,
        'Total_Biallelic_Samples' : 0,
        'Total_Sig_Biallelic_Samples' : 0,
        'Total_Sig_ASE_Ref': 0,
        'Total_Sig_ASE_Alt': 0,
        'Total_Biallelic_No_ASE': 0,
        'Total_Passing_Homozygous': 0,
        'Total_Passing_Homozygous_Ref': 0,
        'Total_Passing_Homozygous_Alt': 0,
        'Total_Non_Testable': 0}  

    # Getting 
    meta_analysis_sig_variants_list = final_meta_analysis_sig_variants_results['meta_analysis_sig_variants_list']

    # creating new filtered file
    sig_meta_file_name = meta_analysis_output_folder + "/sig_meta_analysis_variants.txt"

    # Create output folder for failing data for verification purposes
    not_sig_meta_file_name = meta_analysis_output_folder + "/not_sig_meta_analysis_variants.txt"

    # Open the files for writing
    sig_meta__file = open(sig_meta_file_name,'w')
    not_sig_meta_file = open(not_sig_meta_file_name,'w')
    
    #opening up variant file
    testable_variants_file = open(testable_variants_file_name, 'r')

    # reading file line by line
    for line in testable_variants_file:

    #Getting the header information from the file
        if line.startswith('#CHROM'):
    
            #Print header of file to new data_file
            sig_meta__file.write(line)
            not_sig_meta_file.write(line)
            continue
              
        else:
            #Tallying total testable variants
            meta_global_counts_dict['Biallelic_Testable_Variants']+= 1
            
            #Splitting the data
            parsed_line = line.rstrip().split('\t')

            #Getting Variant Results Only
            variant_name = parsed_line[0] + ":" + parsed_line[1]

            if variant_name in meta_analysis_sig_variants_list:
                sig_meta__file.write(line)
                meta_global_counts_dict['Total_Sig_ASE_Variants']+= 1

            else:
                not_sig_meta_file.write(line)

    # Close the files
    sig_meta__file.close()
    not_sig_meta_file.close()

    return{'sig_meta_file_name': sig_meta_file_name, 'meta_global_counts_dict':meta_global_counts_dict} 


def merge_final_results(meta_analysis_output_folder, sig_variants_final_output_file_name, final_meta_analysis_sig_variants_results):

    """

    Functions sorts through the sig_variants file and merges in the meta results data
    so a user can easily examine p-values, q-values and overall sample relationships

    : Param meta_analysis_output_folder: Output location of meta-data
    : Param sig_variants_final_output_file_name: Name of the significant variants to merge with meta data
    : Param final_meta_analysis_sig_variants_results: Meta-data significant variants results from analysis

    : Return NONE:

    """
    
    # Get the meta-analysis results dictionary from analysis
    meta_variants_results_dict = final_meta_analysis_sig_variants_results['meta_variants_results_dict']
    
    # open up the temporary file for analysis (reads data)
    sig_variants_meta_file = open((meta_analysis_output_folder + '/sig_variants_report_and_meta_results.txt'),'w')
    
    # open up significant variants results file (reads data)
    sig_variants_file = open(sig_variants_final_output_file_name,'r')

    # Parsing through the data
    for line in sig_variants_file:

        #Skipping the header information from the file
        if line.startswith('#CHROM'):
            line = line.rstrip('\n')
            
            sig_variants_meta_file.write(line + "\t\t\tAnalyzed_Values\tChi_Square_Value\t" +
                                         "Degrees_of_Freedom\tMeta_p-value\tFDR_q-value\tSignificant\n")

        #go through the rest of the lines
        else:

            #Removing the new line character
            line = line.rstrip('\n')

            #parse the data line
            parsed_line=line.split('\t')

            #variant name (used for tallying significant variants)
            variant=(parsed_line[0])+":"+(parsed_line[1])

            if variant in meta_variants_results_dict:
                analyzed_values = meta_variants_results_dict[variant]['analyzed_values']
                chi_square_value = meta_variants_results_dict[variant]['chi_square_value']
                degrees_of_freedom = meta_variants_results_dict[variant]['degrees_of_freedom']
                meta_p_value = meta_variants_results_dict[variant]['original_p_value']
                meta_q_value = meta_variants_results_dict[variant]['q_value']
                meta_verdict = meta_variants_results_dict[variant]['verdict']

                sig_variants_meta_file.write(line + "\t\tmeta_results\t" + str(analyzed_values)
                                             + "\t" + str(chi_square_value)
                                             + "\t" + str(degrees_of_freedom)
                                             + "\t" + str(meta_p_value)
                                             + "\t" + str(meta_q_value)
                                             + "\t" + str(meta_verdict) + "\n")
                
            else:
                continue
                
    # Closing the files
    sig_variants_meta_file.close()
    sig_variants_file.close()
    
    return()


def tally_final_meta_results(sig_meta_file_name, meta_global_counts_dict, meta_sample_p_value_cutoff):

    """

    Function parses through and analyzes all the results from the FDR correction
    and tallies the final results for final printing of reports

    : Param sigm_meta_file_name: filename of all significant variants
    : Param meta_global_counts_dict: dictionary of meta-results
    : Param meta_sample_p_value_cutoff: cutoff p-value for tallying samples (estimation)

    : Return samples_list: List of sampels
    : Return samples_counters_dict: Dictionary with information about sample results
    : Return variant_list: List of variants being analyzed
    : Return variant_results_dictionary: Tallying dictionary for variants
    : Return variant_header: Header of file with variants
    : Return tally_global_dictionary: Global tallying dictionary filled in with values

    """

    # Re-name the variable for the function
    tally_global_dictionary = meta_global_counts_dict

    #opening up variant file
    final_results_file = open(sig_meta_file_name, 'r')

    #Create a variant dictionary
    variant_results_dictionary = {}

    #Create a variant list for cycling through
    variant_list = []
    
    #reading file line by line
    for line in final_results_file:

        #Getting the header information from the file
        if line.startswith('#CHROM'):
            
            #Splitting the data
            parsed_line = line.rstrip().split('\t')

            # Getting Variant Header
            variant_header = parsed_line[:8]
            variant_header = ('\t'.join(map(str,variant_header)))

            # Retrieve_Sample_Names_for_Logging
            samples_list = parsed_line[9:]
            numb_individuals = len(samples_list)

            #Create a tallying dictionary for results
            samples_counters_dict = create_sample_tallying_counters(samples_list)

            pass
   
        #Getting Variant Results Only
        else:
            
            #Splitting the data
            parsed_line = line.rstrip().split('\t')
            
            #Variant information
            variant_info = parsed_line[:8]
            variant_key = ('\t'.join(map(str,variant_info)))

            # Add new variant to list
            variant_list.append(variant_key)

            variant_results_dictionary.update({variant_key: {'Biallelic_Testable': 0, 'Sig_ASE': 0, 
                           'Sig_ASE_Ref': 0, 'Sig_ASE_Alt': 0,
                           'Biallelic_No_ASE': 0,
                           'Passing_Homozygous': 0,
                           'Homozygous_Ref': 0,
                           'Homozygous_Alt': 0,                                  
                           'Non_Testable': 0}}) 
            
            # Getting variant results
            variant_sample_results=parsed_line[9:]

            #Starting reading and writing 
            x=0
            for sample_data in variant_sample_results:

                #get sample name from Index Location
                sample = samples_list[x]

                # Tallying all tests possible
                tally_global_dictionary ['Total_Tests'] += 1

                #Splitting the data for later parsing
                split_sample_data = sample_data.split(':')

                if split_sample_data[0]=='Biallelic':
                    tally_global_dictionary ['Total_Biallelic_Samples'] += 1
                    samples_counters_dict[sample]['Biallelic_Testable'] += 1
                    variant_results_dictionary[variant_key]['Biallelic_Testable'] += 1

                    # Getting the p-value of interest from the data (reference index)
                    sample_p_value = float(split_sample_data[3])

                    #Test if breach p_value cutoff
                    if float(sample_p_value) < float(meta_sample_p_value_cutoff):
                        tally_global_dictionary ['Total_Sig_Biallelic_Samples'] += 1
                        samples_counters_dict[sample]['Sig_ASE'] += 1
                        variant_results_dictionary[variant_key]['Sig_ASE'] += 1

                        #Look to see if ref or alt is higher
                        variant_counts = (split_sample_data[2])
                        split_counts = variant_counts.split(',')

                        # If Ref is Higher Than Alt
                        if int(split_counts[0]) > int(split_counts[1]):
                            samples_counters_dict[sample]['Sig_ASE_Ref'] += 1
                            tally_global_dictionary ['Total_Sig_ASE_Ref'] += 1
                            variant_results_dictionary[variant_key]['Sig_ASE_Ref'] += 1

                            # Move Counter
                            x+= 1

                        # If Alt is Higher Than Ref
                        else:
                            samples_counters_dict[sample]['Sig_ASE_Alt'] += 1
                            tally_global_dictionary ['Total_Sig_ASE_Alt'] += 1
                            variant_results_dictionary[variant_key]['Sig_ASE_Alt'] += 1

                            # Move Counter
                            x+= 1
                            
                    #Biallic but not significant (fails cutoff)
                    else:
                        samples_counters_dict[sample]['Biallelic_No_ASE'] += 1
                        tally_global_dictionary ['Total_Biallelic_No_ASE'] += 1
                        variant_results_dictionary[variant_key]['Biallelic_No_ASE'] += 1

                        # Move Counter
                        x+= 1
                        
                        
                #Homozygous passing sample 
                elif split_sample_data[0]=='Homo':
                    samples_counters_dict[sample]['Passing_Homozygous'] += 1
                    tally_global_dictionary ['Total_Passing_Homozygous'] += 1
                    variant_results_dictionary[variant_key]['Passing_Homozygous'] += 1

                    #Tally Type of Homozygous Genotype (allow freq calculations) 
                    homozygous_genotype = (split_sample_data[1])

                    if homozygous_genotype == "0/0":
                        samples_counters_dict[sample]['Passing_Homozygous_Ref'] += 1
                        variant_results_dictionary[variant_key]['Homozygous_Ref'] += 1
                        tally_global_dictionary ['Total_Passing_Homozygous_Ref'] += 1
                        # Move Counter
                        x+= 1

                    elif homozygous_genotype == "1/1":
                        samples_counters_dict[sample]['Passing_Homozygous_Alt'] += 1
                        variant_results_dictionary[variant_key]['Homozygous_Alt'] += 1
                        tally_global_dictionary ['Total_Passing_Homozygous_Alt'] += 1
                        # Move Counter
                        x+= 1

                    else:
                        print ("Danger Incorrect Genotype")
                        print ("Printing line of problematic data")
                        print (line)
                        print ("Printing problematic genotype")
                        print (homozygous_genotype)
                        sys.exit()
                            
                # No passing sample  
                else:
                    samples_counters_dict[sample]['Non_Testable'] += 1
                    tally_global_dictionary['Total_Non_Testable'] += 1
                    variant_results_dictionary[variant_key]['Non_Testable'] += 1

                    # Move Counter
                    x+= 1

    #Closing the initial file
    final_results_file.close()

    return (samples_list, samples_counters_dict, variant_list, variant_results_dictionary, variant_header,
            tally_global_dictionary)


#########################################################################################################################
##################################### Printing Tallying Results  ########################################################
#########################################################################################################################

# Setup generic printing code with standard inputs to allow for this part of code
# to be recycled when performing a different statistical test to try and identify
# ASE among variants


def create_sample_tallying_dict(samples_list):

    """

    Creats a tallying dictionary for the samples used ONLY
    for "testable" variants global file

    : Param samples_list: List of the samples

    : Return tallying_dict: Empty dictionary to store
        all the values

    """

    # Create dictionary
    tallying_dict = {}

    # Loop over the samples list
    for sample in samples_list:

        # Add entry to dictionary
        tallying_dict.update({sample: {'total_count': 0, 'testable_biallelic': 0}})

    return(tallying_dict)


def count_sample_biallelic_testable(input_file_name):


    """

    Counts the number of biallelic testable variants for
    a sample.

    : Param input_file_name: Name of the file being parsed

    : Return testable_samples_counters_dict: Dictionary 


    """

    # Open the file
    input_file = open(input_file_name, 'r')

    # Starting looping over lines of the file
    for line in input_file:

        # Remove the new line for safety
        line = line.rstrip("\n")

        # Getting the header information from the file (sample names)
        if line.startswith('#CHROM'):
            
            # Splitting the data
            parsed_line = line.rstrip().split('\t')

            # Getting Variant Header
            variant_header = parsed_line[:8]
            variant_header = ('\t'.join(map(str,variant_header)))

            # Retrieve_Sample_Names_for_Logging
            samples_list = parsed_line[9:]
            numb_individuals = len(samples_list)

            #Create a tallying dictionary for results
            testable_samples_counters_dict = create_sample_tallying_dict(samples_list)

        # Start Examing the Samples
        else:
            
            # Splitting the data
            parsed_line = line.rstrip().split('\t')

            # Getting Variant Header
            variant_header = parsed_line[:8]
            variant_header = ('\t'.join(map(str,variant_header)))

            # Retrieve sample calls
            sample_allelic_results_calls = parsed_line[9:]

            # Sample index counter
            sample_index_counter = 0

            # Loop over the samples
            for sample_call in sample_allelic_results_calls:

                # Get Sample ID
                sample_id = samples_list[sample_index_counter]

                # Add to the dictionary counter
                testable_samples_counters_dict[sample_id]['total_count'] +=1

                # Get the Biallelic passing variants
                if sample_call.startswith('Biallelic'):
                    testable_samples_counters_dict[sample_id]['testable_biallelic'] +=1

                # Skip over the rest of the sample variant verdicts
                # Possible options: No_Data, Homo_Low_Count, Homo, Low_Read_Count, Low_Allele_Count
                else:
                    pass
                    
                # Add to the counter
                sample_index_counter +=1

    # Close the file
    input_file.close()
   
    return(testable_samples_counters_dict)


def create_sample_tallying_counters(samples_list):
    
    """
    
    Creates a tallyiing dictionary of samples for reporting
    the final results
    
    : Param samples_list: List of samples

    : Return samples_counter_dict: Dictionary of samples and empty tally
        scores
    
    """

    samples_counters_dict = {}
    tallying_dictionary = {'Biallelic_Testable': 0, 'Sig_ASE': 0, 
                           'Sig_ASE_Ref': 0, 'Sig_ASE_Alt': 0,
                           'Biallelic_No_ASE': 0,
                           'Passing_Homozygous': 0,
                           'Passing_Homozygous_Ref': 0,
                           'Passing_Homozygous_Alt': 0,
                           'Non_Testable': 0} 

    for sample in samples_list:

        #Making a deep copy of the dictionary, so each one is independent
        tallying_dictionary = copy.deepcopy(tallying_dictionary)
        
        samples_counters_dict.update({sample: tallying_dictionary})
    
    return(samples_counters_dict)


def printing_sample_results(folder_pathway, samples_list, samples_counters_dict, testable_file_biallelic_samples_dict):
    
    """

    Prints a new file of sample tally reports for ASE
    
    : Param folder_pathway: Pathway of where to put the data after tallying
    : Param samples_list: List of samples to cycle through when printing
    : Param samples_counters_dict: Tallied dictionary of all the ASE results
    
    : Return NONE:
    
    """
    
    samples_report_file = open(folder_pathway + "/sig_samples_report.txt", "w")

    samples_report_file.write("Sample\tTestableFile_Biallelic_Testable\t\tSigResultsFile_Biallelic_Testable\tBiallelic_No_ASE\tSig_ASE\t\tSig_ASE_Ref\tSig_ASE_Alt\t"
                             "\tHomo_Passing\t\tHomo_Ref\tHomo_Alt\t\tNon-Testable\n")

    for sample in samples_list:
        biallelic_testable = str(samples_counters_dict[sample]['Biallelic_Testable'])
        sig_ase = str(samples_counters_dict[sample]['Sig_ASE'])
        sig_ase_ref = str(samples_counters_dict[sample]['Sig_ASE_Ref'])
        sig_ase_alt = str(samples_counters_dict[sample]['Sig_ASE_Alt'])
        biallelic_no_ase = str(samples_counters_dict[sample]['Biallelic_No_ASE'])
        passing_homozygous = str(samples_counters_dict[sample]['Passing_Homozygous'])
        passing_homozygous_ref = str(samples_counters_dict[sample]['Passing_Homozygous_Ref'])
        passing_homozygous_alt = str(samples_counters_dict[sample]['Passing_Homozygous_Alt'])
        non_testable = str(samples_counters_dict[sample]['Non_Testable'])

        # Get the global testable count from the testable file not just from the significant file
        testablefile_biallelic_testable = str(testable_file_biallelic_samples_dict[sample]['testable_biallelic'])
 
        samples_report_file.write(str(sample) + "\t"
                                  + testablefile_biallelic_testable + "\t\t"
                                  + biallelic_testable + "\t"
                                  + biallelic_no_ase + "\t"
                                  + sig_ase + "\t"
                                  + "ASE Breakdown:" + "\t"
                                  + sig_ase_ref + "\t" + sig_ase_alt + "\t\t" 
                                  + passing_homozygous + "\tHomo_Breakdown\t"
                                  + passing_homozygous_ref + "\t" + passing_homozygous_alt
                                  + "\t\t" + non_testable + "\n")
                            
    samples_report_file.close()


def creating_frequency_bin_dict():
    
    '''
    
    Creating a frequency bin to place sample tallies into,
    bins go from 0-9 (example 0 = 0-0.1 frequency)
    Note: If 100%, need to subtract 1 from the value to tally
    correctly in the 0.9-1.0 freq bin.
    
    : Param NONE
    : Return freq_bin_dict: dictionary to place tallying values
    
    '''

    freq_bin_dict ={}
    for x in range(0, 10, 1):

        freq_bin_dict.update({x: 0})

    # Add a final total category for tallying
    freq_bin_dict.update({'total': 0})

    # Tally for number of samples (should be the same for every variant)
    freq_bin_dict.update({'total_sample_count': 0})

    # Making a deep of the dictionary, so each is independent
    freq_bin_dict = copy.deepcopy(freq_bin_dict)

    return(freq_bin_dict)


def printing_frequency_dictionary(folder_pathway, freq_bin_dict, study_file_name,
                                  study_description):
    
    """

    Printing the frequency binning analysis to a file
    
    : Param folder_pathway: Overall location to put the file
    : Param freq_bin_dict: Frequency bin dictionary being printing
    : Param study_file_name: Name of specific binning study
    : Param study_description: Name of study being performed
    
    : Return: NONE
    
    """
    
    dict_out_file = open(folder_pathway + "/" + study_file_name + ".txt", "w")

    dict_out_file.write(study_description + "\n")
    dict_out_file.write("\n")
    dict_out_file.write("Bin\tCounts\n")
    dict_out_file.write("0.0 - 0.1\t" + str(freq_bin_dict[0]) + "\n")
    dict_out_file.write("0.1 - 0.2\t" + str(freq_bin_dict[1]) + "\n")
    dict_out_file.write("0.2 - 0.3\t" + str(freq_bin_dict[2]) + "\n")
    dict_out_file.write("0.3 - 0.4\t" + str(freq_bin_dict[3]) + "\n")
    dict_out_file.write("0.4 - 0.5\t" + str(freq_bin_dict[4]) + "\n")
    dict_out_file.write("0.5 - 0.6\t" + str(freq_bin_dict[5]) + "\n")
    dict_out_file.write("0.6 - 0.7\t" + str(freq_bin_dict[6]) + "\n")
    dict_out_file.write("0.7 - 0.8\t" + str(freq_bin_dict[7]) + "\n")
    dict_out_file.write("0.8 - 0.9\t" + str(freq_bin_dict[8]) + "\n")
    dict_out_file.write("0.9 - 1.0\t" + str(freq_bin_dict[9]) + "\n")
    dict_out_file.write("\n")
    dict_out_file.write("\n")
    dict_out_file.write("Total number of variants analyzed: "
                        + str(freq_bin_dict['total']) + "\n")
    dict_out_file.write("Total number of possible samples per variant: "
                        + str(freq_bin_dict['total_sample_count']))
    dict_out_file.close()
    return()


def printing_variant_results(folder_pathway, variant_list, variant_results_dictionary, variant_header):
    
    """

    Prints a variant report of all variants that show significant ASE results and also triggers
    the printing of all frequency bin tables of results
    
    : Param folder_pathway: Pathway where to output the results
    : Param variant_list: List of variants to keep dictionary in order
    : Param variant_results_dictionary: Dictionary of all tallying of variant data
    : Param variant_header: Varaint header information to output to file
    
    : Return sig_variants_final_output_file_name: Returns file name for meta-analysis
        and gets merged in with the final meta data for sig variants
    
    """

    sig_variants_final_output_file_name = folder_pathway + "/sig_variants_report.txt"
    
    variants_report_file = open(sig_variants_final_output_file_name, "w")

    variants_report_file.write(variant_header+ "\tAlt_Allele_Freq\tBiallelic_Testable\tBiallelic_No_ASE\tSig_ASE\t\tSig_ASE_Ref\tSig_ASE_Alt\t"
                             "\tHomo_Passing\t\tHomo_Ref\tHomo_Alt\t\tNon-Testable\n")

    # Sample Prevalance of ASE Dictionary
    sample_prev_freq_bin_dict = creating_frequency_bin_dict()

    # ASE variant freq bin
    ase_variant_freq_bin_dict = creating_frequency_bin_dict()
    
    for variant in variant_list:
        biallelic_testable = str(variant_results_dictionary[variant]['Biallelic_Testable'])
        sig_ase = str(variant_results_dictionary[variant]['Sig_ASE'])
        sig_ase_ref = str(variant_results_dictionary[variant]['Sig_ASE_Ref'])
        sig_ase_alt = str(variant_results_dictionary[variant]['Sig_ASE_Alt'])
        biallelic_no_ase = str(variant_results_dictionary[variant]['Biallelic_No_ASE'])
        passing_homozygous = str(variant_results_dictionary[variant]['Passing_Homozygous'])

        # Breakdown of Homozygous sample counts 
        homozygous_ref_samples = str(variant_results_dictionary[variant]['Homozygous_Ref'])
        homozygous_alt_samples = str(variant_results_dictionary[variant]['Homozygous_Alt'])

        # Non-testable Sample Counts
        non_testable = str(variant_results_dictionary[variant]['Non_Testable'])

        # Total Sample Count
        total_sample_count = (int(biallelic_testable) + int(passing_homozygous) + int(non_testable))

        ################################################################
        # Analyzing Overall ASE Prevalance Amoung Samples for a Variant
        # Prevlance Amoung Variant
        prevalance = int(int(sig_ase) / total_sample_count * 10)
        # Check not 100%, otherwise subtract 1
        if prevalance == 10:
            prevalance = prevalance - 1
        else:
            pass
        #Adding to frequency bin
        sample_prev_freq_bin_dict[prevalance] += 1
        sample_prev_freq_bin_dict['total'] += 1
        # Just recording for printing purposes (NO TALLYING OCCURRING)
        sample_prev_freq_bin_dict['total_sample_count'] = total_sample_count
        ################################################################

        # Calculating MAF
        # Counts from samples
        biallelic_count = int((variant_results_dictionary[variant]['Biallelic_Testable']))
        homozygous_ref_count = 2*(int(variant_results_dictionary[variant]['Homozygous_Ref']))
        homozygous_alt_count = 2*(int(variant_results_dictionary[variant]['Homozygous_Alt']))

        # Reference Counts
        reference_count = biallelic_count + homozygous_ref_count
        # Alternative Counts
        alternative_count = biallelic_count + homozygous_alt_count
        # Total Counts
        total_count = reference_count + alternative_count                    

        # Getting the minor allele frequency from the counts
        # NOTE: MAF variable name is confusing SHOULD HAVE BEEN- alt_allele_freq
        MAF = str(round((alternative_count / total_count),2))

        ################################################################
        # Analyzing Overall ASE Frequency (What freq shows ASE most)
        # Prevlance Amoung Variant
        freq_bin = int(float(MAF)*10)
        # Check not 100%, otherwise subtract 1
        if freq_bin == 10:
            freq_bin = freq_bin - 1
        else:
            pass
        #Adding to frequency bin
        ase_variant_freq_bin_dict[freq_bin] += 1
        ase_variant_freq_bin_dict['total'] += 1
        # Just recording for printing purposes (NO TALLYING OCCURRING)
        ase_variant_freq_bin_dict['total_sample_count'] = total_sample_count
        ################################################################

        # Printing all variants results to the file                    
        variants_report_file.write(str(variant) + "\t" + MAF + "\t" + biallelic_testable + "\t"
                                  + biallelic_no_ase + "\t"
                                  + sig_ase + "\t"
                                  + "ASE_Breakdown:" + "\t"
                                  + sig_ase_ref + "\t" + sig_ase_alt + "\t\t" 
                                  + passing_homozygous + "\tHomo_Breakdown\t"  + homozygous_ref_samples
                                  + "\t" + homozygous_alt_samples + "\t\t" + non_testable + "\n")

    variants_report_file.close()

    #Printing prevelance frequencying binning results to file
    prev_study_descript = 'Prevalance of an ASE Sample Among ALL Samples for a Variant (freq binning)'
    prevalance_study_file_name = 'ase_freq_prevalence_among_samples'
    print_prevelance_freq = printing_frequency_dictionary(folder_pathway, sample_prev_freq_bin_dict,
                                                          prevalance_study_file_name, prev_study_descript)
    
    #Printing Overall frequencying of a variant identified with ASE
    freq_study_descript = 'Binning of Overall Variant alternative allele frequency which Has At Least One ASE Hit'
    freq_study_file_name = 'alt_allele_freq_binning_one_ase_hit'
    print_sig_var_freq = printing_frequency_dictionary(folder_pathway, ase_variant_freq_bin_dict,
                                                          freq_study_file_name, freq_study_descript)

    return(sig_variants_final_output_file_name)



def printing_summary_report(folder_pathway, tally_global_dictionary,
                            parameter_stuff, statistical_test, raw_rna_seq_stats,
                            avg_ref_allele_ratio, indel_stats_dict,
                            python_version, program_name, date, vadt_output_directory,
                            user_defined_output_location):

    """

    Function prints a summary report file of all the results from the analysis
    process

    : Param folder_pathway: Pathway of output folder (where to put report)
    : Param tally_global_dictionary: Tallying of global counts for summmary report
    : Param parameter_stuff: Original user input paramters passed to the program
    : Param statistical_test: Name of statistical test being used to for statistical correction
    : Param raw_rna_seq_stats: Filtering statistics from the analysis process
    : Param avg_ref_allele_ratio: Average reference allele ratio (calculated from the reference allele bias analysis)
    : Param indel_stats_dict: Dictionary of indel statistics for printing
    : Param python_version: version of python running
    : Param program_name: name of program running
    : Param date: date and time program ran
    : Param vadt_output_directory: name of output directory
    : Param user_defined_output_location: original location defined by the user

    : Return summary_report_file_name: return report file, to record total time of program

    """

    # Create a name of the output report file
    summary_report_file_name = (folder_pathway + "/Summary_Report_" + statistical_test + ".txt")

    # Start Printing Report
    summary_report = open(summary_report_file_name, "w")

    summary_report.write("Log Report of Statistical Analysis Performed\n")
    summary_report.write("\n")
    summary_report.write("\n")
    summary_report.write("Running using Python version: " + python_version + "\n")
    summary_report.write("\n")
    summary_report.write("Program Name is: " + program_name + "\n")
    summary_report.write("\n")
    summary_report.write("Program ran at: " + date + "\n")
    summary_report.write("\n")
    summary_report.write("\n")

    summary_report.write("######################################################################\n")
    summary_report.write("########################## Required Modules ##########################\n")
    summary_report.write("######################################################################\n")
    summary_report.write("\n")
    summary_report.write("from __future__ import division\n")
    summary_report.write("import os.path\n")
    summary_report.write("import sys\n")
    summary_report.write("import time\n")
    summary_report.write("from scipy import stats\n")
    summary_report.write("from scipy.stats import combine_pvalues\n")
    summary_report.write("import numpy as np\n")
    summary_report.write("import copy\n")
    summary_report.write("from datetime import datetime\n")
    summary_report.write("from platform import python_version\n")
    summary_report.write("\n")
    summary_report.write("\n")

    summary_report.write("######################################################################\n")
    summary_report.write("##################### USER PARAMETER SETTINGS ########################\n")
    summary_report.write("######################################################################\n")

    summary_report.write("\n")
    summary_report.write("The input file is: " + str(parameter_stuff['file_name'])+"\n")
    summary_report.write("\n")
    summary_report.write("User defined output location: " + user_defined_output_location + "\n")
    summary_report.write("\n")
    summary_report.write("The minimum qualtity score (phred score) for a variant is: "
                         + str(parameter_stuff['quality_score_min']) + "\n")
    summary_report.write("The indel exclusion region length is: "
                         + str(parameter_stuff['indel_exclusion_region_length']) + " basepairs from an indel\n")
    summary_report.write("The number of allowable reference alleles is (currently program is limited to one): "
                         + str(parameter_stuff['numb_ref_alleles_allowed']) + "\n")
    summary_report.write("The number of allowable alternative alleles is (currently program is limited to one): "
                         + str(parameter_stuff['numb_alt_alleles_allowed']) + "\n")
    summary_report.write("The minimum number of total read counts for a sample per variant is: "
                         + str(parameter_stuff['min_total_read_count']) + "\n")
    summary_report.write("\n")
    summary_report.write("IMPORTANT TO NOTE BINOMIAL PROBABILITY VALUE\n")
    summary_report.write("The binomial probability is set to: "
                         + str(parameter_stuff['binomial_probability_value']) + "\n")
    summary_report.write("\n")
    summary_report.write("\n")

    summary_report.write("######################################################################\n")
    summary_report.write("############## Pathway Info for Output of Program ####################\n")
    summary_report.write("######################################################################\n")

    summary_report.write("\n")
    summary_report.write("User defined output location: " + user_defined_output_location + "\n")
    summary_report.write("\n")
    summary_report.write("Directory created by VADT output: " + vadt_output_directory + "\n")
    summary_report.write("\n")
    summary_report.write("The FULL output directory for analysis is: "
                         + str(parameter_stuff['output_file_location'])+"\n")
    summary_report.write("\n")
    summary_report.write("\n")
    
    summary_report.write("######################################################################\n")
    summary_report.write("###################### Data about file ###############################\n")
    summary_report.write("######################################################################\n")

    summary_report.write("\n")
    summary_report.write("The number of variants in the file: "+str(raw_rna_seq_stats['variants_in_file']) + "\n")
    summary_report.write("The number of samples in the VCF file are: "+str(raw_rna_seq_stats['no_of_individuals_in_file']) + "\n")
    summary_report.write("\n")
    summary_report.write("\n")

    summary_report.write("######################################################################\n")
    summary_report.write("################### Global Filtering of Data for ASE #################\n")
    summary_report.write("######################################################################\n")

    summary_report.write("\n")
    summary_report.write("The number of variants that failed GATK filters: "
                         + str(raw_rna_seq_stats['failed_GATK_SNP_filter']) + "\n")
    summary_report.write("The number of variants that failed Qualtity Score Filter (<20): "
                         + str(raw_rna_seq_stats['failed_Qual_filter']) + "\n")
    summary_report.write("\n")
    summary_report.write("The number of indel exclusion regions identified: "
                         + str(raw_rna_seq_stats['no_indel_exclusion_regions']) + "\n")
    summary_report.write("The number of indels identified (some variants have multiple indels): "
                         + str(indel_stats_dict['number_indels_identified']) + "\n")
    summary_report.write("The longest indel was : "
                         + str(indel_stats_dict['longest_indel']) + " bases \n")
    summary_report.write("The shortest indel was : "
                         + str(indel_stats_dict['shortest_indel']) + " bases \n")
    summary_report.write("The average indel length was : "
                         + str(round(indel_stats_dict['average_indel'], 2)) + " bases \n")
    summary_report.write("The number of variants excluded because of indel exclusion regions: "
                         + str(raw_rna_seq_stats['numb_SNPs_excl_Indels']) + "\n")
    summary_report.write("\n")
    summary_report.write("The number of variants excluded with too many reference alleles: "
                         + str(raw_rna_seq_stats['ref_allele_mutiple_forms']) + "\n")
    summary_report.write("The number of variants excluded with too many alternative alleles: "
                         + str(raw_rna_seq_stats['alt_allele_mutiple_forms']) + "\n")
    summary_report.write("\n")

    summary_report.write("######################################################################\n")
    summary_report.write("############### Sample Level Filtering of Variants ###################\n")
    summary_report.write("######################################################################\n")
    summary_report.write("\n")

    summary_report.write("Each sample is investigated and if ALL samples fail the variant fails\n")
    summary_report.write("\n")
    summary_report.write("The number of variants excluded with no genotype values (all samples): "
                         + str(raw_rna_seq_stats['no_genotype_values']) + "\n")
    summary_report.write("The number of variants excluded with low read counts values (all samples): "
                         + str(raw_rna_seq_stats['low_read_count']) + "\n")
    summary_report.write("The number of variants excluded with low freq count (all samples): "
                         + str(raw_rna_seq_stats['low_freq_count']) + "\n")
    summary_report.write("The number of variants excluded that were all homozygous reference calls (all samples): "
                         + str(raw_rna_seq_stats['no_ref_homozygous_variants']) + "\n")
    summary_report.write("The number of variants excluded that were all homozygous alternative calls (all samples): "
                         + str(raw_rna_seq_stats['no_alt_homozygous_variants']) + "\n")
    summary_report.write("The number of variants excluded that were  a combination of homozygous calls (all samples): "
                         + str(raw_rna_seq_stats['no_combo_homozygous_variants']) + "\n")
    summary_report.write("The number of variants excluded that fail for a combo of filter failures (combination of all prior filters): "
                         + str(raw_rna_seq_stats['combo_filter_failure']) + "\n")
    summary_report.write("\n")
    summary_report.write("The number of passing variants: "
                         + str(raw_rna_seq_stats['passing_variants']) + "\n")
    summary_report.write("\n")
    summary_report.write("The average reference allele ratio for testable variants is: "
                         + str(round(avg_ref_allele_ratio, 4)) + "\n")
    summary_report.write("\n")
    summary_report.write("\n")

    summary_report.write("######################################################################\n")
    summary_report.write("#################### ANALYSIS PERFORMED ##############################\n")
    summary_report.write("######################################################################\n")

    # Getting the p-value cutoffs from the parameter dictionary
    multi_dim_adjust_pvalue_cutoff = parameter_stuff['multi_dim_adjust_pvalue_cutoff']
    meta_BH_adj_p_value_cutoff = parameter_stuff['meta_BH_adj_p_value_cutoff']
    meta_sample_p_value_cutoff = parameter_stuff['meta_sample_p_value_cutoff']

    # Meta-Analysis Report
    if statistical_test == 'meta_analysis':
        summary_report.write("\n")
        summary_report.write("STATISTICAL TEST IMPLEMENTENT: Meta-Analysis with BH P-value Adjustment\n")
        summary_report.write("\n")
        summary_report.write("BH Adjusted P-value Cutoff is: " + str(meta_BH_adj_p_value_cutoff) + "\n")
        summary_report.write("\n")
        summary_report.write("P-value Cutoff for Estimated Sample Tallying: " + str(meta_sample_p_value_cutoff) + "\n")
        summary_report.write("\n")
        summary_report.write("\n")

    # Multi-Dimensional P-value Adjustment
    else:
        summary_report.write("\n")
        summary_report.write("STATISTICAL TEST IMPLEMENTENT: Multi-Dimensional P-Value Adjustment\n")
        summary_report.write("\n")
        summary_report.write("Multi-Dimensional P-value Adjustment Cutoff: " + str(multi_dim_adjust_pvalue_cutoff) + "\n")
        summary_report.write("\n")
        summary_report.write("\n")

    summary_report.write("               SUMMARY REPORT OF RESULTS")
    summary_report.write("\n")
    summary_report.write("The total number of testable biallelic variants is: "        
                         + str(tally_global_dictionary['Biallelic_Testable_Variants']) + "\n")

    summary_report.write("The total number of significant ASE variants is: "
                         + str(tally_global_dictionary['Total_Sig_ASE_Variants']) + "\n")

    summary_report.write("\n")

    summary_report.write("            Breakdown of Significant Variants Results\n")

    summary_report.write("The total number of ALL tests possible: "
                         + str(tally_global_dictionary['Total_Tests']) + "\n")

    summary_report.write("The total number of biallelic variants in all samples: "
                         + str(tally_global_dictionary['Total_Biallelic_Samples']) + "\n")

    summary_report.write("The total number of significant biallelic variants in all samples: "
                         + str(tally_global_dictionary['Total_Sig_Biallelic_Samples']) + "\n")

    summary_report.write("The total number of biallelic variants in all samples with NO significant ASE: "
                         + str(tally_global_dictionary['Total_Biallelic_No_ASE']) + "\n")

    summary_report.write("\n")

    summary_report.write("       Breakdown of Significant Biallelic Variants in All Samples\n")     

    summary_report.write("The total number of significant ASE where reference allele is higher: "
                         + str(tally_global_dictionary['Total_Sig_ASE_Ref']) + "\n")

    summary_report.write("The total number of significant ASE where alternative allele is higher: "
                         + str(tally_global_dictionary['Total_Sig_ASE_Alt']) + "\n")

    summary_report.write("\n")


    summary_report.write("              Breakdown of Homozygous Results\n")

    summary_report.write("The total number of homozygous variants in all samples (pass qualitity controls): "
                         + str(tally_global_dictionary['Total_Passing_Homozygous']) + "\n")

    summary_report.write("The total number of reference homozygous counts: "
                         + str(tally_global_dictionary['Total_Passing_Homozygous_Ref']) + "\n")

    summary_report.write("The total number of alternative homozygous counts: "
                         + str(tally_global_dictionary['Total_Passing_Homozygous_Alt']) + "\n")

    summary_report.write("\n")

    summary_report.write("                  Non-testable Results\n")
    
    summary_report.write("The total number of variants in all samples that were not testable: "
                         + str(tally_global_dictionary['Total_Non_Testable']) + "\n")

    summary_report.write("\n")
    summary_report.write("\n")
    summary_report.write("######################################################################\n")
    summary_report.write("######################################################################\n")
    summary_report.write("######################################################################\n")

    summary_report.close()

    return(summary_report_file_name)


def get_list_sig_variants(file_name):
    
    """

    Parses through the significant variants file
    and creates a dictionary of variants and their
    information for later parsing
    
    : Param file_name: Name of file being parsed

    : Return sig_variants_dict: A dictionary of significant results  

    """

    sig_variants_dict = {}

    sig_variants_file = open(file_name, 'r')

    # Skipping the header
    for line in sig_variants_file:
        if line.startswith('#CHROM'):
            continue
        
        # Getting the actual data to create a dictionary
        else:
            parsed_line = line.rstrip().split('\t')

            variant_name = parsed_line[0] + ":" + parsed_line[1]
            variant_info = parsed_line[:8]
            sig_variants_dict.update({variant_name: variant_info})

    # Closing the file
    sig_variants_file.close()

    return(sig_variants_dict)


def meta_data_for_mapping_bias_plots(testable_variants_file, sig_variants_dict, output_directory):
    
    """

    Function parses through all the testable biallelic variants
    and tallies the overall read counts to investigate reference allele bias
    
    : Param testable_variants_file: File with ALL biallelic testable variants
    : Param sig_variants_dict: All ASE SNPs identified in the study
    : Param output_directory: Tells where to output the data file

    : Return avg_ref_allele_ratio: Average reference allele ratio ---help identify ref allele bias in data
    
    """

    #opening up variant file
    testable_variants_file=open(testable_variants_file, 'r') 

    #Creating a log file for variant meta analysis
    ref_bias_file=open(output_directory + "/data_for_ref_allele_bias_plotting.txt",'w')

    # Writing the heading to file for plotting file
    ref_bias_file.write("Variant_Name\tStatus\tRef_Count\tAlt_Count\tTotal_Count\tRatio\n")

    # Creating Ref Allele Ratio List
    ref_allele_ratio_list = []
    
    #reading file line by line
    for line in testable_variants_file:

        #Getting the header information from the file
        if line.startswith('#CHROM'):
            continue

        else:
            #Splitting the data
            parsed_line = line.rstrip().split('\t')

            variant_name = parsed_line[0] + ":" + parsed_line[1]

            # Determining if the variant showed ASE
            if variant_name in sig_variants_dict:
                status = 'Sig_ASE'
            else:
                status = 'No_ASE'
            

            # Simple Counters
            variant_total_reference_count = 0
            variant_total_alternative_count = 0
            variant_total_count = 0 

            #Getting Variant Results Only
            variant_sample_results=parsed_line[9:]

            #Starting reading and writing 
            for sample_data in variant_sample_results:
                split_sample_data = sample_data.split(':')

                if split_sample_data[0]=='Biallelic':

                    counts = split_sample_data[2].split(',')

                    reference_count = int(counts[0])
                    alternative_count = int(counts[1])
                    total_count = reference_count + alternative_count

                    #add to variant totals
                    variant_total_reference_count += reference_count
                    variant_total_alternative_count += alternative_count
                    variant_total_count += total_count 
 
                else:
                    continue

            # Get the reference allele ratio
            ref_allele_ratio = (variant_total_reference_count/variant_total_count)

            ref_allele_ratio_list.append(ref_allele_ratio)

            # Printing alll results to file
            ref_bias_file.write(str(variant_name) + "\t" + status + "\t" + str(variant_total_reference_count)
                                + "\t" + str(variant_total_alternative_count)
                                + "\t" + str(variant_total_count) + "\t"
                                + str(ref_allele_ratio) + "\n" )

    #Getting the average ref allele ratio
    avg_ref_allele_ratio = np.mean(ref_allele_ratio_list)

    # Closing the file
    testable_variants_file.close()
    ref_bias_file.close()

    return(avg_ref_allele_ratio)
    



#########################################################################################################################
############################################ Main Function ##############################################################
#########################################################################################################################

def main():

    print ("")
    print ("########### VCF ASE Detection Tool (VADT) ###########")
    print ("")
    print ("Starting analysis of data")

    ###############################################################################
    ################ Parsing Through User Parameters ##############################

    # Setting up program
    start_time=time.clock()
    home_directory = os.getcwd()
    
    # Program first test to see if the parameter file is present
    verdict=os.path.isfile('VADT_Parameter_File.txt')
    
    # If Program finds the parameter file it then parses it
    if verdict == True:
        print ("")
        print ("Loading Input from User Provided Parameter File")
        #Code to pull in data from parameter file (TESTING ONLY)
        print ("")
        parameter_stuff = parsing_input_parameter_file('VADT_Parameter_File.txt')

    
    ###Else code pulls in data from the qsub file/command line
    else: 
        user_input=sys.argv[1:]
        user_input=' '.join(user_input)
        print ("")
        print ("Loading input from User Supplied Parameters from Submission File")
        print ("")
        print (user_input)
        parameter_stuff=parsing_input(user_input)
                             
    # Version of python
    python_version = sys.version

    # Program name
    program_name = sys.argv[0]
                             
    # Embed Data and Time (Output Folder)
    date = datetime.now().strftime('%Y-%m-%d' + "_" + '%H-%M')

    # Directory name
    vadt_output_directory = "VADT_output_" + date

    # Making a directory to put all important, but non-essential results in
    # See if directory exists otherwise make it
    output_file_location = parameter_stuff['output_file_location']

    # Keep record of original output location
    user_defined_output_location = output_file_location

    # Update output file location with time stamp
    output_file_location = output_file_location + "/" + vadt_output_directory
    parameter_stuff['output_file_location'] = output_file_location

    ###############################################################################
    ####################### Filtering VCF FIle ####################################
    ###############################################################################
    
    # Making a directory to put all important, but non-essential results in
    # See if directory exists otherwise make it

    # Get output location from dictionary
    output_file_location = parameter_stuff['output_file_location']
    
    verdict = os.path.exists(output_file_location + '/Filtering_Results')
    if str(verdict) == 'False':
        print ("Creating Log Directory of Filtering Results")
        os.makedirs(output_file_location + '/Filtering_Results')
    else:
        print("Log Directory already exists")
        print("")

    print("")
    print ("Filtering the RNA-Seq VCF Results Data Based on User Parameters")
    # Filter RNA_Seq_Data File before begining matching
    filtering_results = filter_RNA_Seq_Data(parameter_stuff)
    filtered_rna_seq_file_name = filtering_results['filtered_rna_seq_file_name']
    raw_rna_seq_stats = filtering_results ['raw_rna_seq_stats']
    indel_stats_dict = filtering_results ['indel_stats_dict']

    # Renaming filtered RNA-Seq File
    testable_variants_file =filtered_rna_seq_file_name

    # Getting a Sample Tally of all testable biallelic counts (individual sample counts)
    testable_file_biallelic_samples_dict = count_sample_biallelic_testable(testable_variants_file)

    ###############################################################################
    ################ Multidimensional P-value Adjustment Code #####################
    ###############################################################################

    print ("Starting Multi-Dimensional p-value Adjustment of Data")
    
    # Creating Output Folder of Results
    verdict = os.path.exists(output_file_location + '/Multi_Dim_Adj_Results')
    if str(verdict) == 'False':
        print ("Creating Multi-Dimensional Adjustment Results folder")
        os.makedirs(output_file_location + '/Multi_Dim_Adj_Results')
    else:
        print("Multi-Dimensional p-value Adjustment folder already exists")
        print("")

    # Name of the pathway for output
    multi_dim_sample_output_dir = output_file_location + '/Multi_Dim_Adj_Results'

    # Getting the p-value cutoffs
    multi_dim_adjust_pvalue_cutoff = parameter_stuff['multi_dim_adjust_pvalue_cutoff']

    print ("The Multi-Dimensional P-Value Cutoff Is:")
    print (multi_dim_adjust_pvalue_cutoff)
    
    # Determine number of pass FDR pvalues after Bonferroni pooling and BH P-value Adjustment
    passing_FDR_pvalue_results = determine_passing_FDR_pvalues(testable_variants_file, multi_dim_adjust_pvalue_cutoff)
    passing_variant_count = passing_FDR_pvalue_results[0]
    total_variants_analyzed = passing_FDR_pvalue_results[1]

    # Analyze the Data and Determine if p-values Significant Based on Multi-Dimensional Correction
    multi_dim_adj_pvalue_file_name = analyze_variants_for_significance_multi_dim_test(testable_variants_file, passing_variant_count,
                                                                                      total_variants_analyzed, multi_dim_sample_output_dir,
                                                                                      multi_dim_adjust_pvalue_cutoff)

    print ("Done Running Multi-Dimensional p-value Adjustment")
    print ("")
    print ("Tallying Results and Creating Reports")
    print ("")
    
    # Determine Number of Significant Variants Based On Signficant Samples
    filter_sig_multi_dim_results = filter_for_multi_dim_sig_samples(multi_dim_sample_output_dir, multi_dim_adj_pvalue_file_name)
    sig_multi_dimensional_file_name = filter_sig_multi_dim_results['sig_multi_dimensional_file_name']
    multi_dim_global_counts_dict = filter_sig_multi_dim_results['multi_dim_global_counts_dict']

    #Investigating Ref Allele Bias in Testable Variants Data
    multi_dim_sig_variants_dict = get_list_sig_variants(sig_multi_dimensional_file_name)

    # Getting Reference Allele Bias
    ########## ADD IN MATPLOT LIB HERE ##############
    multi_dim_avg_ref_allele_ratio = multi_dim_data_for_mapping_bias_plots(testable_variants_file, multi_dim_sig_variants_dict, 
    															 multi_dim_sample_output_dir)
    #################################################

    # Tallying Final Results from All Significant Samples
    tally_final_multi_dim = tally_final_multi_dim_results(sig_multi_dimensional_file_name, multi_dim_global_counts_dict)
    multi_dim_samples_list = tally_final_multi_dim[0]
    multi_dim_samples_counters_dict = tally_final_multi_dim[1]
    multi_dim_variant_list = tally_final_multi_dim[2]
    multi_dim_variant_results_dictionary = tally_final_multi_dim[3]
    multi_dim_variant_header = tally_final_multi_dim[4]
    tally_multi_dim_global_dictionary =  tally_final_multi_dim[5]
    
    #Printing all the final results to various report files
    printing_sample_results(multi_dim_sample_output_dir, multi_dim_samples_list, multi_dim_samples_counters_dict, testable_file_biallelic_samples_dict)

    sig_variants_multi_dim_final_file_name = printing_variant_results(multi_dim_sample_output_dir, multi_dim_variant_list,
                                                                      multi_dim_variant_results_dictionary, multi_dim_variant_header)

    # Define statistical test variable
    statistical_test = 'multi_dimensional_pvalue_adj'

    # Final Summary Report of All the Stats
    summary_report_file_name_multi_dim = printing_summary_report(multi_dim_sample_output_dir, tally_multi_dim_global_dictionary,
                            parameter_stuff, statistical_test, raw_rna_seq_stats,
                            multi_dim_avg_ref_allele_ratio, indel_stats_dict,
                            python_version, program_name, date, vadt_output_directory, user_defined_output_location)

    print ("Done Running Multi-Dimensional Module of VADT")
    print ("")
    print ("")

    ###############################################################################
    ################Meta Analysis of Data for ASE Testing##########################

    print ("Starting Meta-Analysis of ASE Results")

    # Getting the p-value cutoff ############ FIX CUTOFF NAME ##############################################################
    ###############################################################################################################
    meta_BH_adj_p_value_cutoff = parameter_stuff['meta_BH_adj_p_value_cutoff']
    meta_sample_p_value_cutoff = parameter_stuff['meta_sample_p_value_cutoff']

    ###################################################################
    print ("The meta analysis cutoff p value is:")
    
    # Creating a Sample FDR Output Folder of Results
    verdict = os.path.exists(output_file_location + '/Meta_Analysis_Results')
    if str(verdict) == 'False':
        print ("Creating Meta_Analysis_Results folder")
        os.makedirs(output_file_location + '/Meta_Analysis_Results')
    else:
        print("Meta_Analysis_Results folder already exists")
        print("")

    # Full pathway of Meta Analysis Results
    meta_analysis_output_folder = output_file_location + '/Meta_Analysis_Results'

    # Perform meta-analysis of data to create a intermediate log file
    meta_analysis_results = meta_analysis(meta_analysis_output_folder, testable_variants_file)

    # Re-printing log file of results with corrected q-value
    final_meta_analysis_sig_variants_results = meta_printing_q_values(meta_analysis_output_folder, meta_analysis_results, meta_BH_adj_p_value_cutoff)

    # Filtering final data for significant and non-significant variants for tallying
    filtered_meta_results = filter_for_meta_sig_variants(testable_variants_file, meta_analysis_output_folder,
                                                         final_meta_analysis_sig_variants_results)
    sig_meta_file_name = filtered_meta_results['sig_meta_file_name']
    meta_global_counts_dict = filtered_meta_results['meta_global_counts_dict']

    #Investigating Ref Allele Bias in Testable Variants Data
    meta_fdr_sig_variants_dict = get_list_sig_variants(sig_meta_file_name)

    # Getting Reference Allele Bias
    ########## ADD IN MATPLOT LIB HERE ##############
    meta_avg_ref_allele_ratio = meta_data_for_mapping_bias_plots(testable_variants_file, meta_fdr_sig_variants_dict, meta_analysis_output_folder)

    # Tallying All Signficant Results from Meta-Analysis
    tally_final_meta_analysis = tally_final_meta_results(sig_meta_file_name, meta_global_counts_dict, meta_sample_p_value_cutoff)
    meta_samples_list = tally_final_meta_analysis[0]
    meta_samples_counters_dict = tally_final_meta_analysis[1] 
    meta_variant_list = tally_final_meta_analysis[2] 
    meta_variant_results_dictionary = tally_final_meta_analysis[3] 
    meta_variant_header = tally_final_meta_analysis[4] 
    tally_meta_global_dictionary = tally_final_meta_analysis[5] 

    #Printing all the final results to various report files
    printing_sample_results(meta_analysis_output_folder, meta_samples_list, meta_samples_counters_dict, testable_file_biallelic_samples_dict)

    sig_variants_meta_final_file_name = printing_variant_results(meta_analysis_output_folder, meta_variant_list,
                                                                 meta_variant_results_dictionary, meta_variant_header)

    # Define statistical test variable
    statistical_test = 'meta_analysis'

    # Final Summary Report of All the Stats
    summary_report_file_name_meta = printing_summary_report(meta_analysis_output_folder, tally_meta_global_dictionary,
                                                            parameter_stuff, statistical_test, raw_rna_seq_stats,
                                                            meta_avg_ref_allele_ratio, indel_stats_dict,
                                                            python_version, program_name, date, vadt_output_directory, user_defined_output_location)

    # Merge final meta data results into the significant variants results file
    merge_final_results(meta_analysis_output_folder, sig_variants_meta_final_file_name, final_meta_analysis_sig_variants_results)

    print ("Done Running Meta-Analysis of VADT")
    print ("")
    print ("")

    ###############################################################################
    #################### End of Program Timing Stuff #################################


    print ("")
    print ("")
    print ("Program Done Running")

    total_time = time.clock() - start_time
    print ("")
    print ("")
    print("Program ran for a total of " + str(round(total_time, 2))+" seconds")

    # Write the final times to the summary report files
    
    # Multi-Dimensional P-value Adjustment Report File
    sum_multi_dim_file = open(summary_report_file_name_multi_dim, 'a')
    sum_multi_dim_file.write("\n")
    sum_multi_dim_file.write("\n")
    sum_multi_dim_file.write("Program Ran successfully!!!")
    sum_multi_dim_file.write("\n")
    sum_multi_dim_file.write("Total Runtime of program ran for a total of " + str(round(total_time, 2))+" seconds\n")
    sum_multi_dim_file.write("\n")
    sum_multi_dim_file.close()

    # Meta Report File
    sum_meta_file = open(summary_report_file_name_meta, 'a')
    sum_meta_file.write("\n")
    sum_meta_file.write("\n")
    sum_meta_file.write("Program Ran successfully!!!")
    sum_meta_file.write("\n")
    sum_meta_file.write("Total Runtime of program ran for a total of " + str(round(total_time, 2))+" seconds\n")
    sum_meta_file.write("\n")
    sum_meta_file.close()

if __name__ == "__main__":
    
    main()





##################################### Version Control #############################
# VADT_beta_3.0.2.py
# - Fixing MAF label in Meta analysis--- proper term is Alt_Allele_Freq NOT MAF (confusing for either major or minor)
#
# VADT_beta_3.0.1.py
# - Updated sample counts to have a tally of all biallelic testable samples based on the testable variants file
#
# VADT_beta_3.0.0.py 
# - Removed the Sample FDR method and replaced with the Multi_Dimensional p-value adjustment (BIG overhall)
# - Re-arranged a lot of code to reflect this change
# - Created a new output directory with date embedding for the program (prevent re-writing results file)
# - Slightly modified the Summary Report files for each statistical test
# - Fixed overall code notes etc.
#
# VADT_beta_2.0.3
# -Update Meta P-values BH method now
#
#VADT_beta_2.0.2.py Original Version Worked Very Well and Stable













