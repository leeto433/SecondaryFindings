#!/bin/bash

# This line allows the option to only look for variants that are of 
known phase only or variants that are both known and unknown phase
export includephase=unphased  # Options are phased or unphased


# export 
vcf=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/GRCh37/VariantCalls/Nexteraexome37Mb/20180518_Nextera37Mb_VariantCalls/20180518_Nextera37Mb_VariantCalls_Split_Annotated_4/20180518_Nextera37Mb_VariantCalls_Split_ann.vcf.gz 
#Creates a vcf variable that indicates which vcf file to use
export ped=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/Meta/Ped.txt # 
Creates a variable that indicates the file that shows pedigree 
information for the cohort
export disorders=/resource/domains/STUDENT/leeto433/diseases3.txt # 
Creates a variable that indicates the file that has the information on 
the diseases that we are looking for


export 
vcf=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/GRCh37/VariantCalls/AV5UTRs/20180505_AV5UTRs_VariantCalls/20180505_AV5UTRs_VariantCalls_Split_Annotated_1/20180505_AV5UTRs_VariantCalls_Split_ann.vcf.gz
export samples="3075,3447,3690,3412"

#export 
vcf=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/GRCh37/VariantCalls/SeqCapEZ2/20180517_SeqCapEZ2_VariantCalls/20180517_SeqCapEZ2_VariantCalls_Split_Annotated_1/20180517_SeqCapEZ2_VariantCalls_Split_ann.vcf.gz
#export samples="1203"

export project=${SCRATCH}/$(basename ${vcf} .vcf.gz)_SecondaryFindings # 
Creates a variable so that we can create a directory that is located in 
the scratch folder of the person running the script, and removes the 
.vcf.gz suffix and appends _SecondaryFindings onto the end of the 
directory name
mkdir -p ${project} # Creates the directory named as above. 'mkdir -p' 
means make parent directories as needed
cd $project # Go to the directory that was just created
mkdir -p slurm # Create a directory for slurm output

# specifies samples to analyse for secondary findings, set samples to 
empty to analyse all the samples
# export samples="2973" 
if [[ ! -z ${samples ]]; then
	echo ${samples} | tr "," "\n" > ${project}/selectedsamples.list
	inputsamples=${project}/selectedsamples.list
else
	inputsamples=${project}/allsamples.list	
fi


module purge # Clear any modules that may already have been loaded. Some 
modules may interfere with BCFtools
module load BCFtools # BCFtools is a program that allows us to work with 
vcf files

# If the allsamples.list file does not already exist, then create the 
file which contains a list of the sample IDs. The -f flag tests whether 
the file exists and is a regular file
if [ ! -f $project/allsamples.list ]; then
	$(which bcftools) query -l ${vcf} > $project/allsamples.list # 
query -l prints list of sample IDs only
fi

