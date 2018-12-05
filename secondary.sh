#!/bin/bash

BASEDIR=$(dirname $0)

# This line allows the option to only look for variants that are of known phase only or variants that are both known and unknown phase
export includephase=unphased  # Options are phased or unphased


# Creates a variable that indicates the file that shows pedigree information for the cohort
export ped=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/Meta/Ped.txt 

# Creates a variable that indicates the file that has the information on the diseases that we are looking for
export disorders=/resource/domains/STUDENT/leeto433/diseases3.txt 

#Creates a vcf variable that indicates which vcf file to use
export vcf=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/GRCh37/VariantCalls/AV5UTRs/20180505_AV5UTRs_VariantCalls/20180505_AV5UTRs_VariantCalls_Split_Annotated_1/20180505_AV5UTRs_VariantCalls_Split_ann.vcf.gz
export samples="3075,3447,3690,3412"
# export vcf=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/GRCh37/VariantCalls/SeqCapEZ2/20180517_SeqCapEZ2_VariantCalls/20180517_SeqCapEZ2_VariantCalls_Split_Annotated_1/20180517_SeqCapEZ2_VariantCalls_Split_ann.vcf.gz
# export samples="1203"
# export vcf=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/GRCh37/VariantCalls/Nexteraexome37Mb/20180518_Nextera37Mb_VariantCalls/20180518_Nextera37Mb_VariantCalls_Split_Annotated_4/20180518_Nextera37Mb_VariantCalls_Split_ann.vcf.gz 
# export samples="2973" 

# if the email variable has not been previously been specified, then try to find the email address for the the user from the /etc/slurm/userlist.txt file
if [ -z ${email} ]; then
	oldIFS=$IFS
	IFS=$'\n'
	userList=($(cat /etc/slurm/userlist.txt | grep $USER))
	for entry in ${userList[@]}; do
		testUser=$(echo $entry | awk -F':' '{print $1}')
		if [ "$testUser" == "$USER" ]; then
			export email=$(echo $entry | awk -F':' '{print $3}')
			break
		fi
	done
	IFS=$oldIFS
	if [ -z ${email} ]; then
		(echo "FAIL: Unable to locate email address for $USER in /etc/slurm/userlist.txt!" 1>&2)
		exit 1
	else
		export MAIL_TYPE=REQUEUE,FAIL,END
		(printf "%-22s%s (%s)\n" "Email address" "${email}" "$MAIL_TYPE" 1>&2)
	fi
fi
if [[ ${email} != "none" ]]; then
	mailme="--mail-user ${email} --mail-type $MAIL_TYPE"
fi

# Creates a variable so that we can create a directory that is located in the scratch folder of the person running the script, and removes the .vcf.gz suffix and appends _SecondaryFindings onto the end of the directory name
export project=${SCRATCH}/$(basename ${vcf} .vcf.gz)_SecondaryFindings 
# Creates the directory named as above. 'mkdir -p' means make parent directories as needed
mkdir -p ${project} 
# Go to the directory that was just created
cd $project 
# Create a directory for slurm output
mkdir -p slurm 

# Clear any modules that may already have been loaded. Some modules may interfere with BCFtools
# BCFtools is a program that allows us to work with vcf files
module purge 
# BCFtools is a program that allows us to work with vcf files
module load BCFtools
# If the allsamples.list file does not already exist, then create the file which contains a list of the sample IDs. The -f flag tests whether the file exists and is a regular file
if [ ! -f $project/allsamples.list ]; then
	# query -l in bcftools prints list of sample IDs only
	$(which bcftools) query -l ${vcf} > ${project}/allsamples.list 
fi

# specifies samples to analyse for secondary findings, set samples to empty to analyse all the samples
if [[ ! -z ${samples} ]]; then
	for i in $(echo ${samples} | tr "," " "); do
		if ! grep -qx "${i}" ${project}/allsamples.list; then
			echo -e "${i} cannot be found in ${vcf}"
			exit
		fi
	done
	echo ${samples} | tr "," "\n" > ${project}/selectedsamples.list
	inputsamples=${project}/selectedsamples.list
else
	inputsamples=${project}/allsamples.list	
fi

# Clear any modules that may already have been loaded. Some modules may interfere with BCFtools
# BCFtools is a program that allows us to work with vcf files
module purge 
# BCFtools is a program that allows us to work with vcf files
module load BCFtools
# If the allsamples.list file does not already exist, then create the file which contains a list of the sample IDs. The -f flag tests whether the file exists and is a regular file
if [ ! -f $project/allsamples.list ]; then
	# query -l in bcftools prints list of sample IDs only
	$(which bcftools) query -l ${vcf} > ${project}/allsamples.list 
fi
# make an array of the sample IDs that will be processed
export SAMPLEARRAY=($(cat ${inputsamples} | tr "\n" " "))
# the number of entries in samplearray - should be the number of samples to be analysed
export NUMSAMPLES=${#SAMPLEARRAY[@]}

#sbatch -J Secondary_Sample_Analysis ${mailme} --array 1-${NUMSAMPLES}%6 ${BASEDIR}/secondary2.sl

echo "${SAMPLEARRAY[@]}"
echo "${NUMSAMPLES}"
echo "${mailme}"
