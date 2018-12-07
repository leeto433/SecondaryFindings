#!/bin/bash

BASEDIR=$(dirname $0)
defaultdisorders=${BASEDIR}/ACMG2.0.txt


# Prompts user to specify the VCF file
while [ -z ${vcf} ] || [[ ! -f ${vcf} ]]; do 
	echo -e "Specify the full path to a compressed VCF file that contains sequence data for your subjects.\n"
	read -e -p "Provide the full path, or q to quit, and press [RETURN]: " vcf
	if [[ ${vcf} == "q" ]]; then exit; fi
	if [ ! -z ${vcf} ] && [[ ! -f ${vcf} ]]; then 
		echo -e "\n**********\nCan't find ${vcf}\n**********\n"
	fi
	if [ -z ${vcf} ]; then 
		echo -e "\n**********\nNo file specified\n**********\n"
	fi
done
echo -e "\nThe VCF file being used is $(basename ${vcf})"
export vcf



if [ -z ${disorders} ]; then
	disorders=walrus
elif [[ ! -f ${disorders} ]]; then
	echo -e "Can't find ${disorders}"
	exit
fi
# Prompts user to specify file containing diseases and genes to search for
while [ ${disorders} == 'walrus' ] || [[ ! -f ${disorders} ]]; do 
	echo -e "\nThe default disease gene table is $(basename ${defaultdisorders}). If you want to use the default file, then enter nothing."
	read -e -p "If you want to use an alternative file, then specify the full path, or q to quit, and press [RETURN]: " disorders
	if [[ ${disorders} == "q" ]]; then exit; fi
	if [ ! -z ${disorders} ] && [[ ! -f ${disorders} ]]; then 
		echo -e "\n**********\nCan't find ${disorders}\n**********\n"
	fi
	if [ -z ${disorders} ]; then 
		echo -e "\n**********\nNo file specified\nThe default file will be used\n**********\n"
		disorders=${defaultdisorders}
	fi	
done
echo -e "\nThe disease gene table being used is $(basename ${disorders})"
export disorders



if [[ ! -z ${ped} ]] && [[ ! -f ${ped} ]]; then 
	echo -e "Cant't find ${ped}"
	exit
fi
# Prompts user to specify pedigree file
while [ -z ${ped} ] || [[ ! -f ${ped} ]]; do 
	echo -e "\nSpecify the full path to the pedigree file\n"
	read -e -p "Provide the full path, or q to quit, and press [RETURN]: " ped
	if [[ ${ped} == "q" ]]; then exit; fi
	if [ ! -z ${ped} ] && [[ ! -f ${ped} ]]; then 
		echo -e "\n**********\nCan't find ${ped}\n**********\n"
	fi
	if [ -z ${ped} ]; then 
		echo -e "\n**********\nNo file specified\n**********\n"
	fi	
done 
echo -e "\nThe pedigree file being used is $(basename ${ped})"
export ped

# Clear any modules that may already have been loaded. Some modules may interfere with BCFtools
# BCFtools is a program that allows us to work with vcf files
module purge 
# BCFtools is a program that allows us to work with vcf files
module load BCFtools
# Creates a variable that contains sample ID list in vcf file
allsamples=$($(which bcftools) query -l ${vcf})


if [ -z ${samples} ]; then 
	samples=walrus
elif [ -f ${samples} ]; then
	for i in $(cat ${samples}); do 
		if ! $(echo -e "${allsamples}" | grep -qx "${i}"); then
			echo -e "\nThe specified sample ${i} from the file ${samples} cannot be found in ${vcf}"
			exit
		fi
	done
else 
	for i in $(echo "${samples}" | tr "," " "); do 
		if ! $(echo -e "${allsamples}" | grep -qx "${i}"); then
			echo -e "\nThe specified sample ${i} from your list cannot be found in ${vcf}"
			exit
		fi
	done
fi



# Prompts user to specify samples 
while [[ ${samples} == walrus ]]; do 
	read -e -p $'\nTo specify subjects to analyse, provide their IDs seperated by commas. If the sample IDs are in a file, then provide the full path. Otherwise, leave blank to analyse all subjects in the VCF file and press [RETURN]: ' samples
	if [[ ${samples} == "q" ]]; then 
		exit
	elif [ -f ${samples} ] && [ ! -z ${samples} ]; then
		for i in $(cat ${samples}); do 
			if ! $(echo -e "${allsamples}" | grep -qx "${i}"); then
				echo -e "\nThe specified sample ${i} from the file ${samples} cannot be found in ${vcf}"
				samples=walrus
				break
			fi
		done
	
	elif [ ! -z ${samples} ]; then
		for i in $(echo "${samples}" | tr "," " "); do 
			if ! $(echo -e "${allsamples}" | grep -qx "${i}"); then
				echo -e "\nThe specified sample ${i} from your list cannot be found in ${vcf}"
				samples=walrus
				break
			fi
		done
	fi
done




# This line allows the option to only look for variants that are of known phase only or variants that are both known and unknown phase
while [[ ${includephase} != "yes" ]] && [[ ${includephase} != "no" ]]; do 
	echo -e "\nDo you wish to include unphased variants when considering potential compound heterozygotes in recessive disorders?"
	# This may need a bit more explanation
	read -e -p "Type yes to include unphased variants, and type no to exclude unphased variants, or press q to quit and press [RETURN]: " includephase
	if [[ ${includephase} == "q" ]]; then exit; fi
done
export includephase # yes means include unphased. no means do not include unphased variants. 


# Creates a variable that indicates the file that shows pedigree information for the cohort
#ped=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/Meta/Ped.txt 

# Creates a variable that indicates the file that has the information on the diseases that we are looking for
#disorders=/resource/domains/STUDENT/leeto433/diseases3.txt 


#Creates a vcf variable that indicates which vcf file to use
#export vcf=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/GRCh37/VariantCalls/AV5UTRs/20180505_AV5UTRs_VariantCalls/20180505_AV5UTRs_VariantCalls_Split_Annotated_1/20180505_AV5UTRs_VariantCalls_Split_ann.vcf.gz
#export samples="3075,3447,3690,3412"
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
mkdir -p ${project}/slurm 

echo -e "${allsamples}" > ${project}/allsamples.list

if [[ ! -z ${samples} ]] && [[ -f ${samples} ]]; then
	cp ${samples} ${project}/selectedsamples.list
# If they specify subject ID numbers
elif [[ ! -z ${samples} ]]; then
	echo ${samples} | tr "," "\n" > ${project}/selectedsamples.list
else
	echo -e "${allsamples}" > ${project}/selectedsamples.list
fi

inputsamples=${project}/selectedsamples.list
echo -e "The sample IDs being investigated are \n$(cat ${inputsamples})"



# make an array of the sample IDs that will be processed
export SAMPLESTRING=$(cat ${inputsamples} | tr "\n" ",")
SAMPLEARRAY=($(cat ${inputsamples} | tr "\n" " ")) 
# the number of entries in samplearray - should be the number of samples to be analysed
export NUMSAMPLES=${#SAMPLEARRAY[@]}

exit

sbatch -J Secondary_Sample_Analysis ${mailme} --array 1-${NUMSAMPLES}%6 ${BASEDIR}/secondary2.sl

#echo "Sample array is ${SAMPLEARRAY[@]}"
#echo "${NUMSAMPLES}"
#echo "${mailme}"
