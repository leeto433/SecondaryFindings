#!/bin/bash

BASEDIR=$(dirname $0)
defaultdisorders=${BASEDIR}/ACMG2.0.txt

function usage () {
echo -e "\

\nThis script creates slurm jobs to analyse samples for secondary findings. 

usage: $0 options:

-h 	Print full help/usage information.

Optional:

-v	Specify a VCF file.
	Provide the full path. 

-d 	Specify a file containing disease gene table.
	Provide the full path. 
	To use the default disease gene table following ACMG2.0 guidelines, 
	type ACMG.

-p 	Specify a file containing pedigree information.
	Provide the full path.

-s 	Specify individual samples to be analysed by entering sample IDs seperated by commas.
	If sample IDs are in a file, then provide the full path.
	If you wish to analyse all the samples, then type all.  
	
-i 	To include unphased variants when considering compound heterozygotes for recessive diseases,
	type yes. 
	To exclude unphased variants, type no. 
	
-b	To include compound variants using haplotype aware BCFtools BCSQ caller, 
	type yes.
	To exclude them, type no. 
"
}

while getopts ":v:d:p:s:i:hb:" OPTION; do 
	case $OPTION in 
		v) 	
			export vcf="${OPTARG}"
			;;
		d) 
			export disorders="${OPTARG}"
			;;
		p)
			export ped="${OPTARG}"
			;;
		s)
			export samples="${OPTARG}"
			;;
		i)
			export includephase="${OPTARG}"
			;;
		b)
			export haplotype="${OPTARG}"
			;;
		h)
			usage
			exit 0
			;;
		\?)
			echo -e "\n**************\nInvalid option: ${OPTARG}\nFor valid options, use -h option\n**************\n"
			exit 0
			;;
		:)
			echo -e "\n**************\nOption ${OPTARG} requires an argument\n**************\n"
			exit 0
			;;
	esac
done	


# Prompts user to specify the VCF file
<<<<<<< HEAD
while [ -z ${vcf} ] || [[ ! -f ${vcf} ]]; do 
	echo -e "Specify the full path to a compressed VCF file that contains sequence data for your subjects.\nNote: This may take a few seconds\n"
=======
# while the vcf variable is empty, or vcf variable is not a file, or vcf variable does not have vcf.gz extension. 
# Will only enter while loop if vcf variable hasn't been set, or specified properly yet
while [ -z ${vcf} ] || [[ ! -f ${vcf} ]] || [[ ! $(echo ${vcf} | grep '.vcf.gz$') =~ ".vcf.gz" ]]; do 
	echo -e "Specify the full path to a compressed VCF file that contains sequence data for your subjects.\n"
>>>>>>> 538dd784390d6712391858a3a5e9a1caf222d6b6
	read -e -p "Provide the full path, or q to quit, and press [RETURN]: " vcf
	if [[ ${vcf} == "q" ]]; then exit; fi
	# if vcf variable is set but not a filename
	if [ ! -z ${vcf} ] && [[ ! -f ${vcf} ]]; then 
		echo -e "\n**********\nCan't find ${vcf}\n**********\n"
	fi
	#if vcf variable is set, is a filename but is not a vcf.gz file
	if  [ ! -z ${vcf} ] && [[ -f ${vcf} ]] && [[ ! $(echo ${vcf} | grep '.vcf.gz$') =~ ".vcf.gz" ]]; then
		echo -e "\n**********\nFile specified is not right file type\n**********\n"	
	fi
	# if vcf variable is unset
	if [ -z ${vcf} ]; then 
		echo -e "\n**********\nNo file specified\n**********\n"
	fi
done
echo -e "\nThe VCF file being used is $(basename ${vcf})"
export vcf



<<<<<<< HEAD


=======
if [ -z ${disorders} ]; then
	disorders=walrus
elif [ ${disorders} == "ACMG" ]; then
	disorders=${defaultdisorders}
elif [[ ! -f ${disorders} ]]; then
	echo -e "Can't find ${disorders}"
	exit
fi
>>>>>>> 538dd784390d6712391858a3a5e9a1caf222d6b6
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




<<<<<<< HEAD



=======
>>>>>>> 538dd784390d6712391858a3a5e9a1caf222d6b6
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


<<<<<<< HEAD




# BCFtools is a program that allows us to work with vcf files
module purge 
# BCFtools is a program that allows us to work with vcf files
module load BCFtools
# Creates a variable so that we can create a directory that is located in the scratch folder of the person running the script, and removes the .vcf.gz suffix and appends _SecondaryFindings onto the end of the directory name
export project=${SCRATCH}/$(basename ${vcf} .vcf.gz)_SecondaryFindings 
# Creates the directory named as above. 'mkdir -p' means make parent directories as needed
mkdir -p ${project} 
# Go to the directory that was just created
cd $project 




# If the allsamples.list file does not already exist, then create the file which contains a list of the sample IDs. The -f flag tests whether the file exists and is a regular file
if [ ! -f ${project}/allsamples.list ]; then
	# query -l in bcftools prints list of sample IDs only
	$(which bcftools) query -l ${vcf} > ${project}/allsamples.list 
fi








# Prompts user to specify samples 
read -e -p $'\nTo specify subjects to analyse, provide their IDs seperated by commas. If the sample IDs are in a file, then provide the full path. Otherwise, leave blank to analyse all subjects in the VCF file and press [RETURN]: ' samples

# If they specify subject ID numbers
if [[ ! -z ${samples} ]] && [[ ! -f ${samples} ]]; then
	#Checking that the ID is in VCF file
	for i in $(echo ${samples} | tr "," " "); do
		if ! grep -qx "${i}" ${project}/allsamples.list ; then
			echo -e "\n${i} cannot be found in ${vcf}"
			exit
		fi
	done
	echo ${samples} | tr "," "\n" > ${project}/selectedsamples.list
	inputsamples=${project}/selectedsamples.list
		
# if they have specified a file	
elif [[ ! -z ${samples} ]] && [[ -f ${samples} ]]; then
	#Checking that the ID is in VCF file
	for i in $(cat ${samples}); do 
		if ! grep -qx "${i}" ${project}/allsamples.list; then
			echo -e "\n${i} cannot be found in ${vcf}"
			exit
		fi
	done	
	cat ${samples} > ${project}/selectedsamples.list
	inputsamples=${project}/selectedsamples.list
	
# If they did not enter any information	
elif [ -z ${samples} ]; then
	inputsamples=${project}/allsamples.list	
=======
if [ -z ${samples} ]; then 
	samples=walrus
elif [ -f ${samples} ]; then
	for i in $(cat ${samples}); do 
		if ! $(echo -e "${allsamples}" | grep -qx "${i}"); then
			echo -e "\nThe specified sample ${i} from the file ${samples} cannot be found in ${vcf}"
			exit
		fi
	done
elif [ ${samples} != "all" ]; then
	for i in $(echo "${samples}" | tr "," " "); do 
		if ! $(echo -e "${allsamples}" | grep -qx "${i}"); then
			echo -e "\nThe specified sample ${i} from your list cannot be found in ${vcf}"
			exit
		fi
	done
>>>>>>> 538dd784390d6712391858a3a5e9a1caf222d6b6
fi



<<<<<<< HEAD






#phase for recessive disorders
=======
# Prompts user to specify samples 
while [[ ${samples} == walrus ]]; do 
	read -e -p $'\nTo specify subjects to analyse, provide their IDs seperated by commas. \nIf the sample IDs are in a file, then provide the full path, or q to quit. \nOtherwise, leave blank to analyse all subjects in the VCF file and press [RETURN]: ' samples
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


>>>>>>> 538dd784390d6712391858a3a5e9a1caf222d6b6

# This line allows the option to only look for variants that are of known phase only or variants that are both known and unknown phase
while [[ ${includephase} != "yes" ]] && [[ ${includephase} != "no" ]]; do 
 	echo -e "\nWhen considering variants to be considered as potential compound heterozygotes, do you wish to include unphased variants? \nIn this case, significant variants will be reported if two or more are present in the same gene, \nirrespective of whether they come from different parents."
	echo -e "\nDo you wish to include unphased variants when considering potential compound heterozygotes in recessive disorders?"
	# This may need a bit more explanation
	read -e -p "Type yes to include unphased variants, and type no to exclude unphased variants, or press q to quit and press [RETURN]: " includephase
	if [[ ${includephase} == "q" ]]; then exit; fi
done
if [[ ${includephase} == "yes" ]]; then
	echo -e "\nUnphased variants will be included for consideration as compound heterozygotes\n"
elif [[ ${includephase} == "no" ]]; then
	echo -e "\nUnphased variants will not be included for consideration as compound heterozygotes\n"
fi
export includephase # yes means include unphased. no means do not include unphased variants. 

# Ask if they want to use haplotype aware calls (BCFtools CSQ annotation)
while [[ ${haplotype} != "yes" ]] && [[ ${haplotype} != "no" ]]; do
	echo -e "\nDo you want to include consequences for compound variants, as determined by BCFtools haplotype aware caller that provides BCSQ annotations?"
	read -e -p "Type yes to include compound variants, or no to exclude them, or press q to quit and press [RETURN]: " haplotype
	if [[ ${haplotype} == "q" ]]; then exit; fi
done
if [[ ${haplotype} == "yes" ]]; then
	echo -e "\nCompound variants will be included\n"
elif [[ ${haplotype} == "no" ]]; then
	echo -e "\nCompound variants will be excluded\n"
fi
export haplotype 

# Creates a variable that indicates the file that shows pedigree information for the cohort
<<<<<<< HEAD
#export ped=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/Meta/Ped.txt 

# Creates a variable that indicates the file that has the information on the diseases that we are looking for
#export disorders=/resource/domains/STUDENT/leeto433/diseases3.txt 

=======
#ped=/mnt/hcs/WCHP_Clinical_Genetics/SequenceData/Meta/Ped.txt 

# Creates a variable that indicates the file that has the information on the diseases that we are looking for
#disorders=/resource/domains/STUDENT/leeto433/diseases3.txt 
>>>>>>> 538dd784390d6712391858a3a5e9a1caf222d6b6


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
	mailme="--mail-user=${email} --mail-type=$MAIL_TYPE"
fi


echo -e "${allsamples}" > ${project}/allsamples.list

<<<<<<< HEAD
=======
if [[ ! -z ${samples} ]] && [[ -f ${samples} ]]; then
	cp ${samples} ${project}/selectedsamples.list
# If they specify subject ID numbers
elif [[ ! -z ${samples} ]] && [[ ${samples} != "all" ]]; then
	echo ${samples} | tr "," "\n" > ${project}/selectedsamples.list
else
	echo -e "${allsamples}" > ${project}/selectedsamples.list
fi

export inputsamples=${project}/selectedsamples.list
echo -e "The sample IDs being investigated are \n$(cat ${inputsamples})"

>>>>>>> 538dd784390d6712391858a3a5e9a1caf222d6b6


# make an array of the sample IDs that will be processed
export SAMPLESTRING=$(cat ${inputsamples} | tr "\n" ",")
SAMPLEARRAY=($(cat ${inputsamples} | tr "\n" " ")) 
# the number of entries in samplearray - should be the number of samples to be analysed
export NUMSAMPLES=${#SAMPLEARRAY[@]}

#exit

cmd1="sbatch -J Secondary_Sample_Analysis ${mailme} --array 1-${NUMSAMPLES}%6 ${BASEDIR}/secondary2.sl"
secondary2_job=$(eval ${cmd1} | awk '{print $4}')

#exit

cmd2="sbatch -J Report_${project} --dependency=afterok:${secondary2_job} ${mailme} ${BASEDIR}/report.sl"
eval ${cmd2}

<<<<<<< HEAD
echo "Sample array is ${SAMPLEARRAY[@]}"
echo "${NUMSAMPLES}"
echo "${mailme}"
=======
#echo "Sample array is ${SAMPLEARRAY[@]}"
#echo "${NUMSAMPLES}"
#echo "${mailme}"
>>>>>>> 538dd784390d6712391858a3a5e9a1caf222d6b6
