#!/bin/bash
#SBATCH --job-name	Secondary
#SBATCH --time		02:00:00
#SBATCH --mem		1G
#SBATCH --cpus-per-task	1
#SBATCH --error		slurm/subject-%A_%a-%j.out
#SBATCH --output	slurm/subject-%A_%a-%j.out	
# This script is to run across a VCF and find secondary findings for each individual

SAMPLEARRAY=($(echo "${SAMPLESTRING}" | tr "," " "))
subject=${SAMPLEARRAY[$(( $SLURM_ARRAY_TASK_ID - 1 ))]}



# Updates jobname on Slurm so we can see which subject we are processing
scontrol update jobid=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} jobname=Secondary_${subject}	

mkdir -p ${project}/${subject} # makes a folder for each subject
cd ${project}/${subject}

if [ -f ${project}/${subject}/${subject}.done ]; then 
	exit 0
fi

father=$(cat $ped | awk -F"\t" -v OFS="\t" "\$2 ~ /^$subject\$/ {print \$3;exit;}") 
mother=$(cat $ped | awk -F"\t" -v OFS="\t" "\$2 ~ /^$subject\$/ {print \$4;exit;}") 

if [[ ${mother} == "0" ]]; then
	mother=""
elif ! grep -qx "${mother}" ${project}/allsamples.list; then
	mother=""
fi

if [[ ${father} == "0" ]]; then
	father=""
elif ! grep -qx "${father}" ${project}/allsamples.list; then
	father=""
fi

# build a list of subject, mother and father if available 
samples=${subject}

if [ ! -z ${mother} ]; then
	samples=${samples},${mother}
fi

if [ ! -z ${father} ]; then
	samples=${samples},${father}
fi

# Updates jobname on Slurm so we can see which subject we are processing
scontrol update jobid=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} jobname=Secondary_${subject}_subsettingVCF
	
if [ ! -f ${project}/${subject}/subset.vcf.gz.done ]; then
	cmd="$(which bcftools) view -Oz -o ${project}/${subject}/subset.vcf.gz -s ${samples} ${vcf}"
	echo "${cmd}"
	eval ${cmd} || exit 1$? 
	$(which bcftools) index ${project}/${subject}/subset.vcf.gz && touch ${project}/${subject}/subset.vcf.gz.done
fi

while read -u 3 -r diseasegene;do
	
	# ignore header line, by testing if it begins with ACMG
	if $(echo "${diseasegene}" | grep -q "^ACMG"); then
		continue
	fi
	
	# Assigns columns to variables
	category=$(echo -e "${diseasegene}" | awk -F"\t" '{print $1}') # ACMG disease category
	acmg_disease=$(echo -e "${diseasegene}" | awk -F"\t" '{print $2}') # ACMG disease name
	mim_disease=$(echo -e "${diseasegene}" | awk -F"\t" '{print $3}') # MIM disease phenotype number
	mim_phenotype=$(echo -e "${diseasegene}" | awk -F"\t" '{print $4}')
	hgnc_symbol=$(echo -e "${diseasegene}" | awk -F"\t" '{print $5}')
	mim_gene=$(echo -e "${diseasegene}" | awk -F"\t" '{print $6}')
	ncbi_gene=$(echo -e "${diseasegene}" | awk -F"\t" '{print $7}')
	inheritance=$(echo -e "${diseasegene}" | awk -F"\t" '{print $8}')
	variations=$(echo -e "${diseasegene}" | awk -F"\t" '{print $9}')
	transcripts=$(echo -e "${diseasegene}" | awk -F"\t" '{print $10}')
	minAF=$(echo -e "${diseasegene}" | awk -F"\t" '{print $11}')
	
	if [ -f ${project}/${subject}/${category}/${hgnc_symbol}/${hgnc_symbol}.done ]; then 
		continue
	fi
		
	expandedvariations=$(echo ${variations} | sed 's/KP/Known pathogenic/g' | sed 's/EP/Expected pathogenic/g' | sed 's/MODERATE/Moderate impact/g' | sed 's/LOW/Low impact/g' | sed 's/MODIFIER/Modifier/g')
	expandedinheritance=$(echo ${inheritance} | sed 's/AD/Autosomal dominant/g' | sed 's/SD/Semi-dominant/g' | sed 's/AR/Autosomal recessive/g' | sed 's/XR/X-linked recessive/g' | sed 's/XD/X-linked dominant/g')
			
	# Updates jobname on Slurm so we can see which subject we are processing
	scontrol update jobid=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} jobname=Secondary_${subject}_FilterVariants

	
	mkdir -p ${project}/${subject}/${category}/${hgnc_symbol}
	cd ${project}/${subject}/${category}/${hgnc_symbol}
	

	echo -e "${subject}" > ${project}/${subject}/${category}/${hgnc_symbol}/subject.txt
	if [ ! -z ${mother} ]; then
		echo -e "${mother}" > ${project}/${subject}/${category}/${hgnc_symbol}/mother.txt
	fi
	if [ ! -z ${father} ]; then
		echo -e "${father}" > ${project}/${subject}/${category}/${hgnc_symbol}/father.txt
	fi
	
	# selects lines with KP
	if [[ ${variations} =~ "KP" ]]; then	
		known="((clinVar_nonMNP_CLNDISDB ~ \"OMIM:${mim_phenotype}\" & (clinVar_nonMNP_GENEINFO ~ \":${ncbi_gene}$\" | clinVar_nonMNP_GENEINFO ~ \":${ncbi_gene}|\" | clinVar_nonMNP_GENEINFO ~ \":${ncbi_gene}&\") & clinVar_nonMNP_CLNSIG ~ \"Pathogenic/i\") | (clinVar_MNP_CLNDISDB ~ \"OMIM:${mim_phenotype}\" & (clinVar_MNP_GENEINFO ~ \":${ncbi_gene}$\" | clinVar_MNP_GENEINFO ~ \":${ncbi_gene}|\") & clinVar_MNP_CLNSIG ~ \"Pathogenic/i\") & GT[@subject.txt] = \"alt\")"
	else
		known=""
	fi
	
	
	
	# selects lines with EP 
	if [[ ${variations} =~ "EP" ]]; then	
		if [ ! -z $transcripts ]; then
			
			transcriptsearch=""
		
			for i in $(echo $transcripts | tr "," " "); do 
				if [ -z "${transcriptsearch}" ]; then 
					transcriptsearch="CSQ ~ \"|HIGH|[^|]*|[^|]*|[^|]*|${i}|\""
				else 
					transcriptsearch="${transcriptsearch} | CSQ ~ \"|HIGH|[^|]*|[^|]*|[^|]*|${i}|\"" 
				fi 
			done
			transcriptsearch="((${transcriptsearch}) & GT[@subject.txt] = \"alt\")"
		else
			transcriptsearch="(CSQ ~ \"|HIGH|${hgnc_symbol}|\" & GT[@subject.txt] = \"alt\")"
			
		fi
		expectedBCSQ=""
		# Now look for compound variants using BCSQ
		if [ ${haplotype} == "yes" ]; then	
			expectedBCSQ="(INFO/BCSQ[*] ~ \"^[^\*]*frameshift\.*|${hgnc_symbol}|[^|]*|[^|]*|[^|]*|[^|]*|[^|]*+\" | INFO/BCSQ[*] ~ \"[^\*]*splice_acceptor\.*|${hgnc_symbol}|[^|]*|[^|]*|[^|]*|[^|]*|[^|]*+\" | INFO/BCSQ[*] ~ \"[^\*]*splice_donor\.*|${hgnc_symbol}|[^|]*|[^|]*|[^|]*|[^|]*|[^|]*+\" | INFO/BCSQ[*] ~ \"[^\*]*stop_gained\.*|${hgnc_symbol}|[^|]*|[^|]*|[^|]*|[^|]*|[^|]*+\" | INFO/BCSQ[*] ~ \"[^\*]*start_lost\.*|${hgnc_symbol}|[^|]*|[^|]*|[^|]*|[^|]*|[^|]*+\" | INFO/BCSQ[*] ~ \"[^\*]*stop_lost\.*|${hgnc_symbol}|[^|]*|[^|]*|[^|]*|[^|]*|[^|]*+\" | INFO/BCSQ[*] ~ \"[^\*]*transcript_ablation\.*|${hgnc_symbol}|[^|]*|[^|]*|[^|]*|[^|]*|[^|]*+\" | INFO/BCSQ[*] ~ \"[^\*]*transcript_amplification\.*|[${hgnc_symbol}|[^|]*|[^|]*|[^|]*|[^|]*|[^|]*+\") & FORMAT/BCSQ[@subject.txt] > 0"
		fi
		expected="(${transcriptsearch}${expectedBCSQ})"
	else
		expected=""			
	fi
	
	
	# selects lines with MODERATE impact. MODERATE impact changes amino acid sequence
	if [[ ${variations} =~ "MODERATE" ]]; then	# tests if variations to report includes PP for this gene
		if [ ! -z $transcripts ]; then # the -z flag returns true if variable is empty. This line returns true if $transcripts is not empty
			
			transcriptsearch="" # resets $transcriptsearch to empty
		
			for i in $(echo $transcripts | tr "," " "); do
				if [ -z "${transcriptsearch}" ]; then # returns true if $transcriptsearch is empty, ie this is the first time through the for loop
					transcriptsearch="CSQ ~ \"|MODERATE|[^|]*|[^|]*|[^|]*|${i}|\""
				else 
					transcriptsearch="${transcriptsearch} | CSQ ~ \"|MODERATE|[^|]*|[^|]*|[^|]*|${i}|\"" 
				fi 
			done
			transcriptsearch="((${transcriptsearch}) & GT[@subject.txt] = \"alt\")"
		else
			transcriptsearch="(CSQ ~ \"|MODERATE|${hgnc_symbol}|\" & GT[@subject.txt] = \"alt\")"
			
		fi
		moderate="${transcriptsearch}"
	else 
		moderate=""
	fi
		
	# selects lines with LOW impact
	if [[ ${variations} =~ "LOW" ]]; then	# tests if variations to report includes LOW for this gene
		if [ ! -z $transcripts ]; then # the -z flag returns true if variable is empty. This line returns true if $transcripts is not empty
			
			transcriptsearch="" # resets $transcriptsearch to empty
		
			for i in $(echo $transcripts | tr "," " "); do
				if [ -z "${transcriptsearch}" ]; then # returns true if $transcriptsearch is empty, ie this is the first time through the for loop
					transcriptsearch="CSQ ~ \"|LOW|[^|]*|[^|]*|[^|]*|${i}|\""
				else 
					transcriptsearch="${transcriptsearch} | CSQ ~ \"|LOW|[^|]*|[^|]*|[^|]*|${i}|\"" 
				fi 
			done
			transcriptsearch="((${transcriptsearch}) & GT[@subject.txt] = \"alt\")"
		else
			transcriptsearch="(CSQ ~ \"|LOW|${hgnc_symbol}|\" & GT[@subject.txt] = \"alt\")"
			
		fi
		low="${transcriptsearch}"
	else 
	low=""
	fi
		
	# selects lines with MODIFIER impact
	if [[ ${variations} =~ "MODIFIER" ]]; then	# tests if variations to report includes MODIFIER for this gene
		if [ ! -z $transcripts ]; then # the -z flag returns true if variable is empty. This line returns true if $transcripts is not empty
			
			transcriptsearch="" # resets $transcriptsearch to empty
		
			for i in $(echo $transcripts | tr "," " "); do
				if [ -z "${transcriptsearch}" ]; then # returns true if $transcriptsearch is empty, ie this is the first time through the for loop
					transcriptsearch="CSQ ~ \"|MODIFIER|[^|]*|[^|]*|[^|]*|${i}|\""
				else 
					transcriptsearch="${transcriptsearch} | CSQ ~ \"|MODIFIER|[^|]*|[^|]*|[^|]*|${i}|\"" 
				fi 
			done
			transcriptsearch="((${transcriptsearch}) & GT[@subject.txt] = \"alt\")"
		else
			transcriptsearch="(CSQ ~ \"|MODIFIER|${hgnc_symbol}|\" & GT[@subject.txt] = \"alt\")"
			
		fi
		modifier="${transcriptsearch}"
	else 
	modifier=""
		
	fi
	
	# Changes filter expression based on what variations to report are specified
	if [[ ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]]; then
		expression="${known}"
	elif [[ ! ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]]; then
		expression="${expected}"
	elif [[ ! ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]]; then
		expression="${moderate}"
	elif [[ ! ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]]; then
		expression="${low}"				
	elif [[ ! ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="${modifier}"		
	elif [[ ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]];then
		expression="(${known} | ${expected})"	
	elif [[ ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]]; then
		expression="(${known} | ${moderate})"
	elif [[ ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]]; then
		expression="(${known} | ${low})"	
	elif [[ ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="(${known} | ${modifier})"	
	elif [[ ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]]; then
		expression="(${known} | ${expected} | ${moderate})"
	elif [[ ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]]; then
		expression="(${known} | ${expected} | ${low})"
	elif [[ ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="(${known} | ${expected} | ${modifier})"
	elif [[ ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="(${known} | ${expected} | ${moderate}| ${modifier})"
	elif [[ ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]]; then
		expression="(${known} | ${expected} | ${moderate} | ${low})"
	elif [[ ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="(${known} | ${expected} | ${moderate} | ${low} | ${modifier})"
	elif [[ ! ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]]; then
		expression="(${expected} | ${moderate})"
	elif [[ ! ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]]; then
		expression="(${expected} | ${moderate} | ${low})"
	elif [[ ! ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="(${expected} | ${moderate} | ${low} | ${modifier})"
	elif [[ ! ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]]; then
		expression="(${expected} | ${low})"
	elif [[ ! ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="(${expected} | ${low} | ${modifier})"
	elif [[ ! ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="(${expected} | ${modifier})"
	elif [[ ! ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="(${expected} | ${moderate} | ${modifier})"
	elif [[ ! ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]]; then
		expression="(${moderate} | ${low})"
	elif [[ ! ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="(${moderate} | ${low} | ${modifier})"	
	elif [[ ! ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="(${moderate} | ${modifier})"
	elif [[ ! ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="(${low} | ${modifier})"
	elif [[ ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]]; then
		expression="(${known} | ${moderate} | ${low})"
	elif [[ ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="(${known} | ${moderate} | ${low} | ${modifier})"
	elif [[ ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="(${known} | ${low} | ${modifier})"
	elif [[ ${variations} =~ "KP" ]] && [[ ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="(${known} | ${expected} | ${low} | ${modifier})"
	elif [[ ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ${variations} =~ "MODIFIER" ]]; then
		expression="(${known} | ${moderate} | ${modifier})"
	elif [[ ! ${variations} =~ "KP" ]] && [[ ! ${variations} =~ "EP" ]] && [[ ! ${variations} =~ "MODERATE" ]] && [[ ! ${variations} =~ "LOW" ]] && [[ ! ${variations} =~ "MODIFIER" ]]; then
		echo -e "No variants were specified for ${hgnc_symbol}" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
		continue					
	fi	
	
	# Filters for gnomADexomes_AF and gnomADgenomes_AF if specified on disease file
	minAFfilter=""
	if [ ! -z ${minAF} ]; then
		minAFfilter="(gnomADexomes_AF<=\"${minAF}\" | gnomADgenomes_AF<=\"${minAF}\" | (gnomADexomes_AF=\".\" & gnomADgenomes_AF=\".\")) &&"	
	fi		
	
	
				
	# Filters for annotations from clinVar nonMNP file: OMIM phenotype number, and NCBI gene ID, and classification includes 'pathogenic', and subjects genotype must contain alternate allele. Also specifies input and output file
	cmd="$(which bcftools) filter --exclude '${transcriptsearch}' --include '${minAFfilter} ${expectedBCSQ}' ${project}/${subject}/subset.vcf.gz -Ov -o ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.vcf" 
	echo "${cmd}"
	eval ${cmd} || exit 1$?
		
	if [ ${inheritance} == "AD" ] || [ ${inheritance} == "SD" ] || [ ${inheritance} == "XD" ]; then	
		
		# Generating information for the report
		# bcftools query tool is extracting information from the vcf file and writing it to a .txt file
		# The expression 'CSQ=%CSQ=CSQ' is an attempt to define the beginning and end of the CSQ annotation for future processing			
		$(which bcftools) query -f '%CHROM\t%POS\t%REF\t%ALT\t%gnomADexomes_AF\t%gnomADgenomes_AF\t%clinVar_nonMNP_VariantType\t%clinVar_nonMNP_CLNHGVS\t%clinVar_nonMNP_CLNSIG\t%clinVar_MNP_VariantType\t%clinVar_MNP_CLNHGVS\t%clinVar_MNP_CLNSIG\tCSQ=%CSQ=CSQ\tBCSQ=%BCSQ=BCSQ[\t%TGT %AD]\n' ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.vcf > ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.txt

		if [ -s ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.txt ]; then
			echo -e "${acmg_disease} ${hgnc_symbol} ${expandedinheritance} ${expandedvariations}" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
		fi

		while read -u 4 -r line; do 
			
			
		
			# Write a header for this variant
			echo -e "CHROM\tPOS\tREF\tALT\tgnomADexomes AF\tgnomADgenomes AF\tClinVar_nonMNP_VariantType\tClinVar_nonMNP_CLNHGVS\tClinVar_nonMNP_CLNSIG\tClinVar_MNP_VariantType\tClinVar_MNP_CLNHGVS\tClinVar_MNP_CLNSIG\t$($(which bcftools) query -l ${project}/${subject}/subset.vcf.gz | sed 's/^'${subject}'$/'${subject}' (Subject)/' | sed 's/^'${mother}'$/'${mother}' (Mother)/' | sed 's/^'${father}'$/'${father}' (Father)/' | tr '\n' '\t')" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt				
			# This prints the line from the text file to the report file excluding the CSQ and BCSQ annotations
			echo -e "${line}" | sed 's/^\(.*\)\tCSQ=.*=CSQ\t\(.*\)$/\1\t\2/' | sed 's/^\(.*\)\tBCSQ=.*=BCSQ\t\(.*\)$/\1\t\2/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
			#This captures the CSQ annotation as a variable
			csq=$(echo -e "${line}" | sed 's/^.*\tCSQ=\(.*\)=CSQ\tBCSQ=.*=BCSQ\t.*$/\1/')
			# This captures the BCSQ annotation as a variable
			bcsq=$(echo -e "${line}" | sed 's/^.*\tCSQ=.*=CSQ\tBCSQ=\(.*\)=BCSQ\t.*$/\1/')
	
			if [ -z ${transcripts} ]; then
						
				count=0
				# This searches the CSQ annotation for each transcript name and if present, will print selected CSQ information to a new line
				for i in $(echo "${csq}" | tr "," " "); do
					if [ $count == 0 ]; then
						echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
					fi
					echo "${i}" | grep "^[^|]*|[^|]*|[^|]*|${hgnc_symbol}|" | sed 's/^[^|]*|\([^|]*\)|\([^|]*\)|\([^|]*\)|[^|]*|[^|]*|\([^|]*\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5\t\6/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
					let "count++"
				done		
			else
				count=0
				for i in $(echo "${transcripts}" | tr "," " "); do
					if [ $count == 0 ]; then
						echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
					fi
					echo "${csq}" | grep "|${i}|" | sed 's/^.*[^|]*|\([^|]*\)|\([^|]*\)|\([^|]*\)|[^|]*|[^|]*|\('${i}'\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5\t\6/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
					let "count++"
				done		
			fi		
		
			count=0
			# This searches the BCSQ annotation for each transcript name and if present, will print selected BCSQ information to a new line
			for i in $(echo "${bcsq}" | tr "," " "); do
				if [ $count == 0 ]; then
					echo -e "\tHaplotype aware consequences (BCSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
				fi
				echo -e "\t\t${i}" | grep "^[^|]*|${hgnc_symbol}|" | tr '|' '\t' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
				let "count++"
			done	
		done 4< ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.txt
		
		if [ -s ${project}/${subject}/${category}/${hgnc_symbol}/report.txt ]; then
			cat ${project}/${subject}/${category}/${hgnc_symbol}/report.txt >> ${project}/${subject}/report.txt
		fi
	
	
	
	
	
	
	
	
	elif [ ${inheritance} == "XR" ]; then	
		
		# Generating information for the report
		# bcftools query tool is extracting information from the vcf file and writing it to a .txt file
		# The expression 'CSQ=%CSQ=CSQ' is an attempt to define the beginning and end of the CSQ annotation for future processing			
		$(which bcftools) query --include '(GT[@subject.txt] = "AA" | GT[@subject.txt] = "A") && GT ~ "."' -f '%CHROM\t%POS\t%REF\t%ALT\t%gnomADexomes_AF\t%gnomADgenomes_AF\t%clinVar_nonMNP_VariantType\t%clinVar_nonMNP_CLNHGVS\t%clinVar_nonMNP_CLNSIG\t%clinVar_MNP_VariantType\t%clinVar_MNP_CLNHGVS\t%clinVar_MNP_CLNSIG\tCSQ=%CSQ=CSQ\tBCSQ=%BCSQ=BCSQ[\t%TGT %AD]\n' ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.vcf > ${project}/${subject}/${category}/${hgnc_symbol}/homselectedvariants.txt

		
		touch ${project}/${subject}/${category}/${hgnc_symbol}/patselectedvariants.txt
		touch ${project}/${subject}/${category}/${hgnc_symbol}/matselectedvariants.txt
		touch ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt
		
		# Finds the variants that are clearly from the father or the mother and writes them to their respective files
		if [ ! -z ${mother} ] && [ ! -z ${father} ]; then
			$(which bcftools) query --include 'GT[@subject.txt] = "RA" && GT[@mother.txt] = "RA" && GT[@father.txt] = "R"' -f '%CHROM\t%POS\t%REF\t%ALT\t%gnomADexomes_AF\t%gnomADgenomes_AF\t%clinVar_nonMNP_VariantType\t%clinVar_nonMNP_CLNHGVS\t%clinVar_nonMNP_CLNSIG\t%clinVar_MNP_VariantType\t%clinVar_MNP_CLNHGVS\t%clinVar_MNP_CLNSIG\tCSQ=%CSQ=CSQ\tBCSQ=%BCSQ=BCSQ[\t%TGT %AD]\n' ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.vcf > ${project}/${subject}/${category}/${hgnc_symbol}/matselectedvariants.txt
			$(which bcftools) query --include 'GT[@subject.txt] = "RA" && GT[@father.txt] = "A" && GT[@mother.txt] = "RR"' -f '%CHROM\t%POS\t%REF\t%ALT\t%gnomADexomes_AF\t%gnomADgenomes_AF\t%clinVar_nonMNP_VariantType\t%clinVar_nonMNP_CLNHGVS\t%clinVar_nonMNP_CLNSIG\t%clinVar_MNP_VariantType\t%clinVar_MNP_CLNHGVS\t%clinVar_MNP_CLNSIG\tCSQ=%CSQ=CSQ\tBCSQ=%BCSQ=BCSQ[\t%TGT %AD]\n' ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.vcf > ${project}/${subject}/${category}/${hgnc_symbol}/patselectedvariants.txt
			# When both parents are heterozygous, we can't tell which parent the variant came from so it gets written to the 'unknown' file
			$(which bcftools) query --include 'GT[@subject.txt] = "RA" && GT[@father.txt] = "A" && GT[@mother.txt] = "RA"' -f '%CHROM\t%POS\t%REF\t%ALT\t%gnomADexomes_AF\t%gnomADgenomes_AF\t%clinVar_nonMNP_VariantType\t%clinVar_nonMNP_CLNHGVS\t%clinVar_nonMNP_CLNSIG\t%clinVar_MNP_VariantType\t%clinVar_MNP_CLNHGVS\t%clinVar_MNP_CLNSIG\tCSQ=%CSQ=CSQ\tBCSQ=%BCSQ=BCSQ[\t%TGT %AD]\n' ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.vcf > ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt			
			$(which bcftools) query --include 'GT[@subject.txt] = "RA" && (GT[@father.txt] ~ "\." || GT[@mother.txt] ~ "\.")' -f '%CHROM\t%POS\t%REF\t%ALT\t%gnomADexomes_AF\t%gnomADgenomes_AF\t%clinVar_nonMNP_VariantType\t%clinVar_nonMNP_CLNHGVS\t%clinVar_nonMNP_CLNSIG\t%clinVar_MNP_VariantType\t%clinVar_MNP_CLNHGVS\t%clinVar_MNP_CLNSIG\tCSQ=%CSQ=CSQ\tBCSQ=%BCSQ=BCSQ[\t%TGT %AD]\n' ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.vcf > ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt			
		
		else 
			# If we don't have information on parents, write to 'unknown' file
			$(which bcftools) query --include 'GT[@subject.txt] = "RA" && GT ~ "."' -f '%CHROM\t%POS\t%REF\t%ALT\t%gnomADexomes_AF\t%gnomADgenomes_AF\t%clinVar_nonMNP_VariantType\t%clinVar_nonMNP_CLNHGVS\t%clinVar_nonMNP_CLNSIG\t%clinVar_MNP_VariantType\t%clinVar_MNP_CLNHGVS\t%clinVar_MNP_CLNSIG\tCSQ=%CSQ=CSQ\tBCSQ=%BCSQ=BCSQ[\t%TGT %AD]\n' ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.vcf > ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt
		fi
		
		# Making a report for homozygous/hemizygous alternate alleles
		if [ -s ${project}/${subject}/${category}/${hgnc_symbol}/homselectedvariants.txt ]; then
			echo -e "Homozygous/hemizygous variants" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
			while read -u 4 -r line; do 
				# Write a header for this variant
				echo -e "CHROM\tPOS\tREF\tALT\tgnomADexomes AF\tgnomADgenomes AF\tClinVar_nonMNP_VariantType\tClinVar_nonMNP_CLNHGVS\tClinVar_nonMNP_CLNSIG\tClinVar_MNP_VariantType\tClinVar_MNP_CLNHGVS\tClinVar_MNP_CLNSIG\t$($(which bcftools) query -l ${project}/${subject}/subset.vcf.gz | sed 's/^'${subject}'$/'${subject}' (Subject)/' | sed 's/^'${mother}'$/'${mother}' (Mother)/' | sed 's/^'${father}'$/'${father}' (Father)/' | tr '\n' '\t')" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt				

				# This prints the line from the text file to the report file excluding the CSQ and BCSQ annotations
				echo -e "${line}" | sed 's/^\(.*\)\tCSQ=.*=CSQ\t\(.*\)$/\1\t\2/' | sed 's/^\(.*\)\tBCSQ=.*=BCSQ\t\(.*\)$/\1\t\2/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
				#This captures the CSQ annotation as a variable
				csq=$(echo -e "${line}" | sed 's/^.*\tCSQ=\(.*\)=CSQ\tBCSQ=.*=BCSQ\t.*$/\1/')
				# This captures the BCSQ annotation as a variable
				bcsq=$(echo -e "${line}" | sed 's/^.*\tCSQ=.*=CSQ\tBCSQ=\(.*\)=BCSQ\t.*$/\1/')
	
				if [ -z ${transcripts} ]; then
						
					count=0
					# This searches the CSQ annotation for each transcript name and if present, will print selected CSQ information to a new line
					for i in $(echo "${csq}" | tr "," " "); do
						if [ $count == 0 ]; then
							echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						fi
						echo "${i}" | grep "^[^|]*|[^|]*|[^|]*|${hgnc_symbol}|" | sed 's/^[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\([^|]*\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						let "count++"
					done		
				else
					count=0
					for i in $(echo "${transcripts}" | tr "," " "); do
						if [ $count == 0 ]; then
							echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						fi
						echo "${csq}" | grep "|${i}|" | sed 's/^.*[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\('${i}'\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						let "count++"
					done		
				fi		
		
				count=0
				# This searches the BCSQ annotation for each transcript name and if present, will print selected BCSQ information to a new line
				for i in $(echo "${bcsq}" | tr "," " "); do
					if [ $count == 0 ]; then
						echo -e "\tHaplotype aware consequences (BCSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
					fi
					echo -e "\t\t${i}" | grep "^[^|]*|${hgnc_symbol}|" | tr '|' '\t' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
					let "count++"
				done	
			done 4< ${project}/${subject}/${category}/${hgnc_symbol}/homselectedvariants.txt
		fi
		if [ ${includephase} == 'yes' ] || [[ -s ${project}/${subject}/${category}/${hgnc_symbol}/matselectedvariants.txt  &&  -s ${project}/${subject}/${category}/${hgnc_symbol}/patselectedvariants.txt ]]; then
		
			# Writing maternal variants to report file
			if [ -s ${project}/${subject}/${category}/${hgnc_symbol}/matselectedvariants.txt ] && [[ $(($(wc -l < ${project}/${subject}/${category}/${hgnc_symbol}/patselectedvariants.txt) + $(wc -l < ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt))) > 0 ]]; then
				echo -e "Maternal variants" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
				while read -u 4 -r line; do 
					# Write a header for this variant
					echo -e "CHROM\tPOS\tREF\tALT\tgnomADexomes AF\tgnomADgenomes AF\tClinVar_nonMNP_VariantType\tClinVar_nonMNP_CLNHGVS\tClinVar_nonMNP_CLNSIG\tClinVar_MNP_VariantType\tClinVar_MNP_CLNHGVS\tClinVar_MNP_CLNSIG\t$($(which bcftools) query -l ${project}/${subject}/subset.vcf.gz | sed 's/^'${subject}'$/'${subject}' (Subject)/' | sed 's/^'${mother}'$/'${mother}' (Mother)/' | sed 's/^'${father}'$/'${father}' (Father)/' | tr '\n' '\t')" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt				

					# This prints the line from the text file to the report file excluding the CSQ and BCSQ annotations
					echo -e "${line}" | sed 's/^\(.*\)\tCSQ=.*=CSQ\t\(.*\)$/\1\t\2/' | sed 's/^\(.*\)\tBCSQ=.*=BCSQ\t\(.*\)$/\1\t\2/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
					#This captures the CSQ annotation as a variable
					csq=$(echo -e "${line}" | sed 's/^.*\tCSQ=\(.*\)=CSQ\tBCSQ=.*=BCSQ\t.*$/\1/')
					# This captures the BCSQ annotation as a variable
					bcsq=$(echo -e "${line}" | sed 's/^.*\tCSQ=.*=CSQ\tBCSQ=\(.*\)=BCSQ\t.*$/\1/')
	
					if [ -z ${transcripts} ]; then
						
						count=0
						# This searches the CSQ annotation for each transcript name and if present, will print selected CSQ information to a new line
						for i in $(echo "${csq}" | tr "," " "); do
							if [ $count == 0 ]; then
								echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							fi
							echo "${i}" | grep "^[^|]*|[^|]*|[^|]*|${hgnc_symbol}|" | sed 's/^[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\([^|]*\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							let "count++"
						done		
					else
						count=0
						for i in $(echo "${transcripts}" | tr "," " "); do
							if [ $count == 0 ]; then
								echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							fi
							echo "${csq}" | grep "|${i}|" | sed 's/^.*[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\('${i}'\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							let "count++"
						done		
					fi		
		
					count=0
					# This searches the BCSQ annotation for each transcript name and if present, will print selected BCSQ information to a new line
					for i in $(echo "${bcsq}" | tr "," " "); do
						if [ $count == 0 ]; then
							echo -e "\tHaplotype aware consequences (BCSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						fi
						echo -e "\t\t${i}" | grep "^[^|]*|${hgnc_symbol}|" | tr '|' '\t' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						let "count++"
					done	
				done 4< ${project}/${subject}/${category}/${hgnc_symbol}/matselectedvariants.txt
			fi
	
			# Writing paternal variants to report file
			if [ -s ${project}/${subject}/${category}/${hgnc_symbol}/patselectedvariants.txt ] && [[ $(($(wc -l < ${project}/${subject}/${category}/${hgnc_symbol}/matselectedvariants.txt) + $(wc -l < ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt))) > 0 ]]; then
				echo -e "Paternal variants" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
				while read -u 4 -r line; do 
					# Write a header for this variant
					echo -e "CHROM\tPOS\tREF\tALT\tgnomADexomes AF\tgnomADgenomes AF\tClinVar_nonMNP_VariantType\tClinVar_nonMNP_CLNHGVS\tClinVar_nonMNP_CLNSIG\tClinVar_MNP_VariantType\tClinVar_MNP_CLNHGVS\tClinVar_MNP_CLNSIG\t$($(which bcftools) query -l ${project}/${subject}/subset.vcf.gz | sed 's/^'${subject}'$/'${subject}' (Subject)/' | sed 's/^'${mother}'$/'${mother}' (Mother)/' | sed 's/^'${father}'$/'${father}' (Father)/' | tr '\n' '\t')" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt				
			
					# This prints the line from the text file to the report file excluding the CSQ and BCSQ annotations
					echo -e "${line}" | sed 's/^\(.*\)\tCSQ=.*=CSQ\t\(.*\)$/\1\t\2/' | sed 's/^\(.*\)\tBCSQ=.*=BCSQ\t\(.*\)$/\1\t\2/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
					#This captures the CSQ annotation as a variable
					csq=$(echo -e "${line}" | sed 's/^.*\tCSQ=\(.*\)=CSQ\tBCSQ=.*=BCSQ\t.*$/\1/')
					# This captures the BCSQ annotation as a variable
					bcsq=$(echo -e "${line}" | sed 's/^.*\tCSQ=.*=CSQ\tBCSQ=\(.*\)=BCSQ\t.*$/\1/')
	
					if [ -z ${transcripts} ]; then
						
						count=0
						# This searches the CSQ annotation for each transcript name and if present, will print selected CSQ information to a new line
						for i in $(echo "${csq}" | tr "," " "); do
							if [ $count == 0 ]; then
								echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							fi
							echo "${i}" | grep "^[^|]*|[^|]*|[^|]*|${hgnc_symbol}|" | sed 's/^[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\([^|]*\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							let "count++"
						done		
					else
						count=0
						for i in $(echo "${transcripts}" | tr "," " "); do
							if [ $count == 0 ]; then
								echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							fi
							echo "${csq}" | grep "|${i}|" | sed 's/^.*[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\('${i}'\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							let "count++"
						done		
					fi		
		
					count=0
					# This searches the BCSQ annotation for each transcript name and if present, will print selected BCSQ information to a new line
					for i in $(echo "${bcsq}" | tr "," " "); do
						if [ $count == 0 ]; then
							echo -e "\tHaplotype aware consequences (BCSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						fi
						echo -e "\t\t${i}" | grep "^[^|]*|${hgnc_symbol}|" | tr '|' '\t' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						let "count++"
					done	
				done 4< ${project}/${subject}/${category}/${hgnc_symbol}/patselectedvariants.txt
			fi
		
			# Writes unknown phase variants to report file
			if [ -s ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt ] && [[ $(($(wc -l < ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt) + $(wc -l < ${project}/${subject}/${category}/${hgnc_symbol}/patselectedvariants.txt) + $(wc -l < ${project}/${subject}/${category}/${hgnc_symbol}/matselectedvariants.txt))) > 1 ]]; then
				echo -e "Variants of unknown phase" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
				while read -u 4 -r line; do 
					# Write a header for this variant
					echo -e "CHROM\tPOS\tREF\tALT\tgnomADexomes AF\tgnomADgenomes AF\tClinVar_nonMNP_VariantType\tClinVar_nonMNP_CLNHGVS\tClinVar_nonMNP_CLNSIG\tClinVar_MNP_VariantType\tClinVar_MNP_CLNHGVS\tClinVar_MNP_CLNSIG\t$($(which bcftools) query -l ${project}/${subject}/subset.vcf.gz | sed 's/^'${subject}'$/'${subject}' (Subject)/' | sed 's/^'${mother}'$/'${mother}' (Mother)/' | sed 's/^'${father}'$/'${father}' (Father)/' | tr '\n' '\t')" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt				
			
					# This prints the line from the text file to the report file excluding the CSQ and BCSQ annotations
					echo -e "${line}" | sed 's/^\(.*\)\tCSQ=.*=CSQ\t\(.*\)$/\1\t\2/' | sed 's/^\(.*\)\tBCSQ=.*=BCSQ\t\(.*\)$/\1\t\2/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
					#This captures the CSQ annotation as a variable
					csq=$(echo -e "${line}" | sed 's/^.*\tCSQ=\(.*\)=CSQ\tBCSQ=.*=BCSQ\t.*$/\1/')
					# This captures the BCSQ annotation as a variable
					bcsq=$(echo -e "${line}" | sed 's/^.*\tCSQ=.*=CSQ\tBCSQ=\(.*\)=BCSQ\t.*$/\1/')
	
					if [ -z ${transcripts} ]; then
						
						count=0
						# This searches the CSQ annotation for each transcript name and if present, will print selected CSQ information to a new line
						for i in $(echo "${csq}" | tr "," " "); do
							if [ $count == 0 ]; then
								echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							fi
							echo "${i}" | grep "^[^|]*|[^|]*|[^|]*|${hgnc_symbol}|" | sed 's/^[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\([^|]*\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							let "count++"
						done		
					else
						count=0
						for i in $(echo "${transcripts}" | tr "," " "); do
							if [ $count == 0 ]; then
								echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							fi
							echo "${csq}" | grep "|${i}|" | sed 's/^.*[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\('${i}'\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							let "count++"
						done		
					fi		
		
					count=0
					# This searches the BCSQ annotation for each transcript name and if present, will print selected BCSQ information to a new line
					for i in $(echo "${bcsq}" | tr "," " "); do
						if [ $count == 0 ]; then
							echo -e "\tHaplotype aware consequences (BCSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						fi
						echo -e "\t\t${i}" | grep "^[^|]*|${hgnc_symbol}|" | tr '|' '\t' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						let "count++"
					done	
				done 4< ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt
			fi
		fi
		
		# Adds heading to report file by inserting disease name, symbol and inheritance to the start of the file
		if [ -s ${project}/${subject}/${category}/${hgnc_symbol}/report.txt ]; then
			sed -i '1s/^/'"${acmg_disease}"' '"${hgnc_symbol}"' '"${expandedinheritance}"' '"${expandedvariations}"'\n/' ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
		fi
		if [ -s ${project}/${subject}/${category}/${hgnc_symbol}/report.txt ]; then
			cat ${project}/${subject}/${category}/${hgnc_symbol}/report.txt >> ${project}/${subject}/report.txt
		fi
	
	
	
	
	
	# If the inheritance is autosomal recessive	
	elif [ ${inheritance} == "AR" ]; then
		# Finds the homozygotes
		$(which bcftools) query --include 'GT[@subject.txt] = "AA" && GT ~ "."' -f '%CHROM\t%POS\t%REF\t%ALT\t%gnomADexomes_AF\t%gnomADgenomes_AF\t%clinVar_nonMNP_VariantType\t%clinVar_nonMNP_CLNHGVS\t%clinVar_nonMNP_CLNSIG\t%clinVar_MNP_VariantType\t%clinVar_MNP_CLNHGVS\t%clinVar_MNP_CLNSIG\tCSQ=%CSQ=CSQ\tBCSQ=%BCSQ=BCSQ[\t%TGT %AD]\n' ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.vcf > ${project}/${subject}/${category}/${hgnc_symbol}/homselectedvariants.txt
		
		touch ${project}/${subject}/${category}/${hgnc_symbol}/patselectedvariants.txt
		touch ${project}/${subject}/${category}/${hgnc_symbol}/matselectedvariants.txt
		touch ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt
		
		
		# Finds the variants that are clearly from the father or the mother and writes them to their respective files
		if [ ! -z ${mother} ] && [ ! -z ${father} ]; then
			cmd="$(which bcftools) query --include 'GT[@subject.txt] = \"RA\" && GT[@mother.txt] = \"RA\" && GT[@father.txt] = \"RR\"' -f '%CHROM\t%POS\t%REF\t%ALT\t%gnomADexomes_AF\t%gnomADgenomes_AF\t%clinVar_nonMNP_VariantType\t%clinVar_nonMNP_CLNHGVS\t%clinVar_nonMNP_CLNSIG\t%clinVar_MNP_VariantType\t%clinVar_MNP_CLNHGVS\t%clinVar_MNP_CLNSIG\tCSQ=%CSQ=CSQ\tBCSQ=%BCSQ=BCSQ[\t%TGT %AD]\n' ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.vcf > ${project}/${subject}/${category}/${hgnc_symbol}/matselectedvariants.txt"
			echo "${cmd}"
			eval ${cmd} || exit 1$?
			$(which bcftools) query --include 'GT[@subject.txt] = "RA" && GT[@father.txt] = "RA" && GT[@mother.txt] = "RR"' -f '%CHROM\t%POS\t%REF\t%ALT\t%gnomADexomes_AF\t%gnomADgenomes_AF\t%clinVar_nonMNP_VariantType\t%clinVar_nonMNP_CLNHGVS\t%clinVar_nonMNP_CLNSIG\t%clinVar_MNP_VariantType\t%clinVar_MNP_CLNHGVS\t%clinVar_MNP_CLNSIG\tCSQ=%CSQ=CSQ\tBCSQ=%BCSQ=BCSQ[\t%TGT %AD]\n' ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.vcf > ${project}/${subject}/${category}/${hgnc_symbol}/patselectedvariants.txt
			
			# When both parents are heterozygous, we can't tell which parent the variant came from so it gets written to the 'unknown' file
			$(which bcftools) query --include 'GT[@subject.txt] = "RA" && GT[@father.txt] = "RA" && GT[@mother.txt] = "RA"' -f '%CHROM\t%POS\t%REF\t%ALT\t%gnomADexomes_AF\t%gnomADgenomes_AF\t%clinVar_nonMNP_VariantType\t%clinVar_nonMNP_CLNHGVS\t%clinVar_nonMNP_CLNSIG\t%clinVar_MNP_VariantType\t%clinVar_MNP_CLNHGVS\t%clinVar_MNP_CLNSIG\tCSQ=%CSQ=CSQ\tBCSQ=%BCSQ=BCSQ[\t%TGT %AD]\n' ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.vcf > ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt			
			$(which bcftools) query --include 'GT[@subject.txt] = "RA" && (GT[@father.txt] ~ "\." || GT[@mother.txt] ~ "\.")' -f '%CHROM\t%POS\t%REF\t%ALT\t%gnomADexomes_AF\t%gnomADgenomes_AF\t%clinVar_nonMNP_VariantType\t%clinVar_nonMNP_CLNHGVS\t%clinVar_nonMNP_CLNSIG\t%clinVar_MNP_VariantType\t%clinVar_MNP_CLNHGVS\t%clinVar_MNP_CLNSIG\tCSQ=%CSQ=CSQ\tBCSQ=%BCSQ=BCSQ[\t%TGT %AD]\n' ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.vcf > ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt			
			# query subject RA, father \. mother RR. Write to paternal 
			# query subject RA, father RR, mother \. Write to maternal
			
									
		# elif mother is present but not the father. Include subject RA, mother RR and write to paternal file. Vice versa for father present but not mother. 
		
		
		
		
		else 
			# If we don't have information on parents, write to 'unknown' file
			$(which bcftools) query --include 'GT[@subject.txt] = "RA" && GT ~ "."' -f '%CHROM\t%POS\t%REF\t%ALT\t%gnomADexomes_AF\t%gnomADgenomes_AF\t%clinVar_nonMNP_VariantType\t%clinVar_nonMNP_CLNHGVS\t%clinVar_nonMNP_CLNSIG\t%clinVar_MNP_VariantType\t%clinVar_MNP_CLNHGVS\t%clinVar_MNP_CLNSIG\tCSQ=%CSQ=CSQ\tBCSQ=%BCSQ=BCSQ[\t%TGT %AD]\n' ${project}/${subject}/${category}/${hgnc_symbol}/selectedvariants.vcf > ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt
		fi
		
		
		
		
		#
		
		
		
		if [ -s ${project}/${subject}/${category}/${hgnc_symbol}/homselectedvariants.txt ]; then
			echo -e "Homozygous variants" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
			while read -u 4 -r line; do 
				# Write a header for this variant
				echo -e "CHROM\tPOS\tREF\tALT\tgnomADexomes AF\tgnomADgenomes AF\tClinVar_nonMNP_VariantType\tClinVar_nonMNP_CLNHGVS\tClinVar_nonMNP_CLNSIG\tClinVar_MNP_VariantType\tClinVar_MNP_CLNHGVS\tClinVar_MNP_CLNSIG\t$($(which bcftools) query -l ${project}/${subject}/subset.vcf.gz | sed 's/^'${subject}'$/'${subject}' (Subject)/' | sed 's/^'${mother}'$/'${mother}' (Mother)/' | sed 's/^'${father}'$/'${father}' (Father)/' | tr '\n' '\t')" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt				

				# This prints the line from the text file to the report file excluding the CSQ and BCSQ annotations
				echo -e "${line}" | sed 's/^\(.*\)\tCSQ=.*=CSQ\t\(.*\)$/\1\t\2/' | sed 's/^\(.*\)\tBCSQ=.*=BCSQ\t\(.*\)$/\1\t\2/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
				#This captures the CSQ annotation as a variable
				csq=$(echo -e "${line}" | sed 's/^.*\tCSQ=\(.*\)=CSQ\tBCSQ=.*=BCSQ\t.*$/\1/')
				# This captures the BCSQ annotation as a variable
				bcsq=$(echo -e "${line}" | sed 's/^.*\tCSQ=.*=CSQ\tBCSQ=\(.*\)=BCSQ\t.*$/\1/')
	
				if [ -z ${transcripts} ]; then
						
					count=0
					# This searches the CSQ annotation for each transcript name and if present, will print selected CSQ information to a new line
					for i in $(echo "${csq}" | tr "," " "); do
						if [ $count == 0 ]; then
							echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						fi
						echo "${i}" | grep "^[^|]*|[^|]*|[^|]*|${hgnc_symbol}|" | sed 's/^[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\([^|]*\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						let "count++"
					done		
				else
					count=0
					for i in $(echo "${transcripts}" | tr "," " "); do
						if [ $count == 0 ]; then
							echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						fi
						echo "${csq}" | grep "|${i}|" | sed 's/^.*[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\('${i}'\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						let "count++"
					done		
				fi		
		
				count=0
				# This searches the BCSQ annotation for each transcript name and if present, will print selected BCSQ information to a new line
				for i in $(echo "${bcsq}" | tr "," " "); do
					if [ $count -eq 0 ]; then
						echo -e "\tHaplotype aware consequences (BCSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
					fi
					echo -e "\t\t${i}" | grep "^[^|]*|${hgnc_symbol}|" | tr '|' '\t' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
					let "count++"
				done	
			done 4< ${project}/${subject}/${category}/${hgnc_symbol}/homselectedvariants.txt
		fi
		
		if [ ${includephase} == 'yes' ] || [[ -s ${project}/${subject}/${category}/${hgnc_symbol}/matselectedvariants.txt  &&  -s ${project}/${subject}/${category}/${hgnc_symbol}/patselectedvariants.txt ]]; then
			if [ -s ${project}/${subject}/${category}/${hgnc_symbol}/matselectedvariants.txt ] && [[ $(($(wc -l < ${project}/${subject}/${category}/${hgnc_symbol}/patselectedvariants.txt) + $(wc -l < ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt))) > 0 ]]; then
				echo -e "Maternal variants" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
				while read -u 4 -r line; do 
					# Write a header for this variant
					echo -e "CHROM\tPOS\tREF\tALT\tgnomADexomes AF\tgnomADgenomes AF\tClinVar_nonMNP_VariantType\tClinVar_nonMNP_CLNHGVS\tClinVar_nonMNP_CLNSIG\tClinVar_MNP_VariantType\tClinVar_MNP_CLNHGVS\tClinVar_MNP_CLNSIG\t$($(which bcftools) query -l ${project}/${subject}/subset.vcf.gz | sed 's/^'${subject}'$/'${subject}' (Subject)/' | sed 's/^'${mother}'$/'${mother}' (Mother)/' | sed 's/^'${father}'$/'${father}' (Father)/' | tr '\n' '\t')" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt				
		
					# This prints the line from the text file to the report file excluding the CSQ and BCSQ annotations
					echo -e "${line}" | sed 's/^\(.*\)\tCSQ=.*=CSQ\t\(.*\)$/\1\t\2/' | sed 's/^\(.*\)\tBCSQ=.*=BCSQ\t\(.*\)$/\1\t\2/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
					#This captures the CSQ annotation as a variable
					csq=$(echo -e "${line}" | sed 's/^.*\tCSQ=\(.*\)=CSQ\tBCSQ=.*=BCSQ\t.*$/\1/')
					# This captures the BCSQ annotation as a variable
					bcsq=$(echo -e "${line}" | sed 's/^.*\tCSQ=.*=CSQ\tBCSQ=\(.*\)=BCSQ\t.*$/\1/')
	
					if [ -z ${transcripts} ]; then
						
						count=0
						# This searches the CSQ annotation for each transcript name and if present, will print selected CSQ information to a new line
						for i in $(echo "${csq}" | tr "," " "); do
							if [ $count -eq 0 ]; then
								echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							fi
							echo "${i}" | grep "^[^|]*|[^|]*|[^|]*|${hgnc_symbol}|" | sed 's/^[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\([^|]*\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							let "count++"
						done		
					else
						count=0
						for i in $(echo "${transcripts}" | tr "," " "); do
							if [ $count -eq 0 ]; then
								echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							fi
							echo "${csq}" | grep "|${i}|" | sed 's/^.*[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\('${i}'\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							let "count++"
						done		
					fi		
		
					count=0
					# This searches the BCSQ annotation for each transcript name and if present, will print selected BCSQ information to a new line
					for i in $(echo "${bcsq}" | tr "," " "); do
						if [ $count -eq 0 ]; then
							echo -e "\tHaplotype aware consequences (BCSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						fi
						echo -e "\t\t${i}" | grep "^[^|]*|${hgnc_symbol}|" | tr '|' '\t' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						let "count++"
					done	
				done 4< ${project}/${subject}/${category}/${hgnc_symbol}/matselectedvariants.txt
			fi
			
		

			if [ -s ${project}/${subject}/${category}/${hgnc_symbol}/patselectedvariants.txt ] && [[ $(($(wc -l < ${project}/${subject}/${category}/${hgnc_symbol}/matselectedvariants.txt) + $(wc -l < ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt))) > 0 ]]; then
				echo -e "Paternal variants" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
				while read -u 4 -r line; do 
					# Write a header for this variant
					echo -e "CHROM\tPOS\tREF\tALT\tgnomADexomes AF\tgnomADgenomes AF\tClinVar_nonMNP_VariantType\tClinVar_nonMNP_CLNHGVS\tClinVar_nonMNP_CLNSIG\tClinVar_MNP_VariantType\tClinVar_MNP_CLNHGVS\tClinVar_MNP_CLNSIG\t$($(which bcftools) query -l ${project}/${subject}/subset.vcf.gz | sed 's/^'${subject}'$/'${subject}' (Subject)/' | sed 's/^'${mother}'$/'${mother}' (Mother)/' | sed 's/^'${father}'$/'${father}' (Father)/' | tr '\n' '\t')" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt				

					# This prints the line from the text file to the report file excluding the CSQ and BCSQ annotations
					echo -e "${line}" | sed 's/^\(.*\)\tCSQ=.*=CSQ\t\(.*\)$/\1\t\2/' | sed 's/^\(.*\)\tBCSQ=.*=BCSQ\t\(.*\)$/\1\t\2/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
					#This captures the CSQ annotation as a variable
					csq=$(echo -e "${line}" | sed 's/^.*\tCSQ=\(.*\)=CSQ\tBCSQ=.*=BCSQ\t.*$/\1/')
					# This captures the BCSQ annotation as a variable
					bcsq=$(echo -e "${line}" | sed 's/^.*\tCSQ=.*=CSQ\tBCSQ=\(.*\)=BCSQ\t.*$/\1/')
	
					if [ -z ${transcripts} ]; then
						
						count=0
						# This searches the CSQ annotation for each transcript name and if present, will print selected CSQ information to a new line
						for i in $(echo "${csq}" | tr "," " "); do
							if [ $count -eq 0 ]; then
								echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							fi
							echo "${i}" | grep "^[^|]*|[^|]*|[^|]*|${hgnc_symbol}|" | sed 's/^[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\([^|]*\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							let "count++"
						done		
					else
						count=0
						for i in $(echo "${transcripts}" | tr "," " "); do
							if [ $count -eq 0 ]; then
								echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							fi
							echo "${csq}" | grep "|${i}|" | sed 's/^.*[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\('${i}'\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							let "count++"
						done		
					fi		
		
					count=0
					# This searches the BCSQ annotation for each transcript name and if present, will print selected BCSQ information to a new line
					for i in $(echo "${bcsq}" | tr "," " "); do
						if [ $count -eq 0 ]; then
							echo -e "\tHaplotype aware consequences (BCSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						fi
						echo -e "\t\t${i}" | grep "^[^|]*|${hgnc_symbol}|" | tr '|' '\t' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						let "count++"
					done	
				done 4< ${project}/${subject}/${category}/${hgnc_symbol}/patselectedvariants.txt
			fi
		
			if [ -s ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt ] && [[  $(($(wc -l < ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt) + $(wc -l < ${project}/${subject}/${category}/${hgnc_symbol}/patselectedvariants.txt) + $(wc -l < ${project}/${subject}/${category}/${hgnc_symbol}/matselectedvariants.txt))) > 1 ]]; then
				echo -e "Variants of unknown phase" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
				while read -u 4 -r line; do 
					# Write a header for this variant
					echo -e "CHROM\tPOS\tREF\tALT\tgnomADexomes AF\tgnomADgenomes AF\tClinVar_nonMNP_VariantType\tClinVar_nonMNP_CLNHGVS\tClinVar_nonMNP_CLNSIG\tClinVar_MNP_VariantType\tClinVar_MNP_CLNHGVS\tClinVar_MNP_CLNSIG\t$($(which bcftools) query -l ${project}/${subject}/subset.vcf.gz | sed 's/^'${subject}'$/'${subject}' (Subject)/' | sed 's/^'${mother}'$/'${mother}' (Mother)/' | sed 's/^'${father}'$/'${father}' (Father)/' | tr '\n' '\t')" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt				
		
					# This prints the line from the text file to the report file excluding the CSQ and BCSQ annotations
					echo -e "${line}" | sed 's/^\(.*\)\tCSQ=.*=CSQ\t\(.*\)$/\1\t\2/' | sed 's/^\(.*\)\tBCSQ=.*=BCSQ\t\(.*\)$/\1\t\2/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
					#This captures the CSQ annotation as a variable
					csq=$(echo -e "${line}" | sed 's/^.*\tCSQ=\(.*\)=CSQ\tBCSQ=.*=BCSQ\t.*$/\1/')
					# This captures the BCSQ annotation as a variable
					bcsq=$(echo -e "${line}" | sed 's/^.*\tCSQ=.*=CSQ\tBCSQ=\(.*\)=BCSQ\t.*$/\1/')
	
					if [ -z ${transcripts} ]; then
						
						count=0
						# This searches the CSQ annotation for each transcript name and if present, will print selected CSQ information to a new line
						for i in $(echo "${csq}" | tr "," " "); do
							if [ $count -eq 0 ]; then
								echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							fi
							echo "${i}" | grep "^[^|]*|[^|]*|[^|]*|${hgnc_symbol}|" | sed 's/^[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\([^|]*\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							let "count++"
						done		
					else
						count=0
						for i in $(echo "${transcripts}" | tr "," " "); do
							if [ $count -eq 0 ]; then
								echo -e "\tVEP Consequences (CSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							fi
							echo "${csq}" | grep "|${i}|" | sed 's/^.*[^|]*|\([^|]*\)|[^|]*|\([^|]*\)|[^|]*|[^|]*|\('${i}'\)|[^|]*|[^|]*|[^|]*|\([^|]*\)|\([^|]*\)|.*$/\t\t\1\t\2\t\3\t\4\t\5/' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
							let "count++"
						done		
					fi		
		
					count=0
					# This searches the BCSQ annotation for each transcript name and if present, will print selected BCSQ information to a new line
					for i in $(echo "${bcsq}" | tr "," " "); do
						if [ $count -eq 0 ]; then
							echo -e "\tHaplotype aware consequences (BCSQ)" >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						fi
						echo -e "\t\t${i}" | grep "^[^|]*|${hgnc_symbol}|" | tr '|' '\t' >> ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
						let "count++"
					done	
				done 4< ${project}/${subject}/${category}/${hgnc_symbol}/unknownselectedvariants.txt
			fi
		fi
		
		# Adds heading to report file by inserting disease name, symbol and inheritance to the start of the file
		if [ -s ${project}/${subject}/${category}/${hgnc_symbol}/report.txt ]; then
			sed -i '1s/^/'"${acmg_disease}"' '"${hgnc_symbol}"' '"${expandedinheritance}"' '"${expandedvariations}"'\n/' ${project}/${subject}/${category}/${hgnc_symbol}/report.txt
		fi
		if [ -s ${project}/${subject}/${category}/${hgnc_symbol}/report.txt ]; then
			cat ${project}/${subject}/${category}/${hgnc_symbol}/report.txt >> ${project}/${subject}/report.txt
		fi
	fi
	touch ${project}/${subject}/${category}/${hgnc_symbol}/${hgnc_symbol}.done	
done 3< ${disorders}
touch ${project}/${subject}/${subject}.done || exit 1

exit 0
# Writes a report to the project level
#if [ -s ${project}/${subject}/report.txt ]; then
#	echo "${subject}" >> ${project}/report.txt
#	cat ${project}/${subject}/report.txt >> ${project}/report.txt
#	echo -e "\n" >> ${project}/report.txt
#fi




