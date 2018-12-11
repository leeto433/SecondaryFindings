#!/bin/bash
#SBATCH --job-name	Secondary_report
#SBATCH --time		00:10:00
#SBATCH --mem		256M
#SBATCH --cpus-per-task	1
#SBATCH --error		slurm/report-%A_%a-%j.out
#SBATCH --output	slurm/report-%A_%a-%j.out	
# This script generates a report at the project level by collating all of the reports from individual samples


if [ ! -f ${project}/${project}_report.txt ]; then 
	echo -e "This is a report of the significant variant analysis of the following individuals: \n\n$(cat ${inputsamples} | tr "\n" "\t")\n\nThe VCF file used was ${vcf}.\nThe disease gene table used was ${disorders}.\nThe pedigree file used was ${ped}." > ${project}/$(basename ${project})_report.txt
	if [ ${includephase} == "yes" ]; then 
		echo -e "CAUTION: This report includes unphased variants for consideration as potential compound heterozygotes.\n" >> ${project}/$(basename ${project})_report.txt
	elif [ ${includephase} == "no" ]; then 
		echo -e "This report does not include unphased variants for consideration as potential compound heterozygotes.\n" >> ${project}/$(basename ${project})_report.txt
	fi
fi

for i in $(cat ${inputsamples} | tr "\n" " "); do
	if [ -s ${project}/${i}/report.txt ]; then
		echo -e "Significant variants for sample ${i}" >> ${project}/$(basename ${project})_report.txt
		cat ${project}/${i}/report.txt >> ${project}/$(basename ${project})_report.txt
		echo -e "\n" >> ${project}/$(basename ${project})_report.txt
	elif [ ! -s ${project}/${i}/report.txt ]; then
		echo -e "Sample ${i} does not have any significant variants\n" >> ${project}/$(basename ${project})_report.txt
	fi
done
exit 0
