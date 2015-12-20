#!/bin/bash

#echo "\\documentclass[a4paPercent_,10pt]{article}" > SimulateGenome.all.stat.tex
#echo "\\usepackage{graphicx}" >> SimulateGenome.all.stat.tex
#echo "\\usepackage{longtable}" >> SimulateGenome.all.stat.tex
#echo "\\begin{document}" >> SimulateGenome.all.stat.tex
#echo "\\includegraphics[scale=0.5]{Rplot_SimulateGenome.pdf}" >> SimulateGenome.all.stat.tex 
#echo "\\newline" >> SimulateGenome.all.stat.tex
#echo "\\begin{tabular}{c|ccc|cc}" >> SimulateGenome.all.stat.tex
#echo "&& Match && \%Mapped& \\\\">>SimulateGenome.all.stat.tex
###### For Genomic #######
#echo "\\begin{longtable}{c|c|c|ccc|c}" >> SimulateGenome.all.stat.tex
#echo "&&Number&&Aligned&& Not Aligned\\\\" >> SimulateGenome.all.stat.tex
#echo "Length& Program & of Reads & Right-position & & Wrong-position\\\\ \hline">> SimulateGenome.all.stat.tex
###### For Bacteria ######
#echo "\\begin{longtable}{c|c|c|c|c}" >> SimulateGenome.all.stat.tex
#echo "&&Number&& Not\\\\" >> SimulateGenome.all.stat.tex
#echo "Length& Program & of Reads & Aligned& Aligned \\\\ \hline">> SimulateGenome.all.stat.tex
#echo "\\endhead" >> SimulateGenome.all.stat.tex

#rm SimulateGenome.all.Rate
#rm "fakeGenome.all.stat"
#Generate Reads
#./readSim -i "fakeGenome.fa" -o "SimulatedGenome" -b 20 -e 50 -g 10000000 -c "fakeGenome" -d 0.9 -u 0.92 -m 0.002 -s 0.8 -v 1 &&

#Indexed Genome (BWA & rCandy)
#bwa index BigChromosomes_whole_genome.fa &&
#/home/public/usr/bin/anfo-index BigChromosomes_whole_genome.fa -o BigChromosomes_whole_genome.anfoIndex &&
#bwa index whole_genome.fa && 
#/home/bioinf/usr/bin/anfo-index  BigChromosomes_whole_genome.fa  -o BigChromosomes_whole_genome.anfoIndex &&

#Do this for all of them
#for length in 25 30
#for length in `seq 25 5 40`
#do
#	for namePrefix in GenomicRead #Bacteria 
#	do
#		echo "Length ${length}"
		#Convert to bam
#		/usr/bin/samtools view /mnt/scratch/Homa/ChoppedReads/${namePrefix}_"$length".sam -Sb > /mnt/scratch/Homa/ChoppedReads/${namePrefix}_"$length".bam
#	done
#done

#for length in `seq 25 5 40`
for length in 25  
do

	for namePrefix in GenomicRead # Bacteria
#	for namePrefix in `splited_*`
	do
		echo "${namePrefix} Length ${length}" 
		#Mapped by BWA ancient parameters
		#time bwa bam2bam -t 40 -g /r1/people/homa_amini/Repository/SimulateGenome/BigChromosomes_whole_genome.fa -n 0.01 -o 2 -f /mnt/scratch/Homa/ChoppedReads/aligned_BWA_WMD_${namePrefix}_${length}.bam /mnt/scratch/Homa/ChoppedReads/${namePrefix}_${length}.bam 1>/mnt/scratch/Homa/ChoppedReads/aligned_BWA_WMD_${namePrefix}_${length}.log 2>/mnt/scratch/Homa/ChoppedReads/aligned_BWA_WMD_${namePrefix}_${length}.err &
		#bwa bam2bam -t 40 -g /r1/people/homa_amini/Repository/SimulateGenome/SimulateGenome.tmp.fa -n 0.01 -o 2 -f /mnt/scratch/Homa/ChoppedReads/aligned_BWA_WMD_${namePrefix}_${length}.bam /mnt/scratch/Homa/ChoppedReads/${namePrefix}_${length}.bam && 
#		/home/homa_amini/Repository/Hash-rCandy/Hash_rCandy -i "/mnt/scratch/Homa/ChoppedReads/aligned_BWA_${namePrefix}_${length}.bam" -o "/mnt/scratch/Homa/ChoppedReads/BWA_All_AS_" -m "BWA-${namePrefix}-0" -l ${length} -r 1 >> SimulateGenome.all.Rate

		#Mapped by BWA default parameters
		#time bwa bam2bam -t 40 -g /r1/people/homa_amini/Repository/SimulateGenome/BigChromosomes_whole_genome.fa  -f /mnt/scratch/Homa/ChoppedReads/aligned_BWA_${namePrefix}_${length}.bam /mnt/scratch/Homa/ChoppedReads/${namePrefix}_${length}.bam 1>/mnt/scratch/Homa/ChoppedReads/aligned_BWA_${namePrefix}_${length}.log 2>/mnt/scratch/Homa/ChoppedReads/aligned_BWA_${namePrefix}_${length}.err &
		#bwa bam2bam -t 40 -g /r1/people/homa_amini/Repository/SimulateGenome/SimulateGenome.tmp.fa  -f /mnt/scratch/Homa/ChoppedReads/aligned_BWA_${namePrefix}_${length}.bam /mnt/scratch/Homa/ChoppedReads/${namePrefix}_${length}.bam &&


		for AS in 20
#		for AS in `seq 10 0.5 20`
		do
			#Mapped with r-candy with parameters
			~udo_stenzel/.cabal/bin/r-candy -o /mnt/scratch/Homa/ChoppedReads/aligned_rCandy_${namePrefix}_${length}_${AS}.bam -A ${AS} -l 0.0  -r 0.0  -d 0.00 -s 0.0 -G /r1/people/homa_amini/Repository/SimulateGenome/BigChromosomes_whole_genome.anfoIndex /mnt/scratch/Homa/ChoppedReads/${namePrefix}_${length}.bam   
			#~/.cabal/bin/r-candy -o /mnt/scratch/Homa/ChoppedReads/aligned_rCandy_${namePrefix}_${length}_${AS}.bam -A ${AS} -l 0.3  -r 0.3  -d 0.02 -s 0.9 -G /r1/people/homa_amini/Repository/SimulateGenome/Genome.Simulated.tmp.anfoIndex /mnt/scratch/Homa/ChoppedReads/${namePrefix}_${length}.bam &  
			#time ~/.cabal/bin/r-candy -o /mnt/scratch/Homa/ChoppedReads/aligned_rCandy_${namePrefix}_${length}.bam -A ${AS} -l 0.3  -r 0.3  -d 0.02 -s 0.9 -G /r1/people/homa_amini/Repository/SimulateGenome/Genome.Simulated.tmp.anfoIndex  /mnt/scratch/Homa/ChoppedReads/${namePrefix}_${length}.bam 1>/mnt/scratch/Homa/ChoppedReads/aligned_rCandy_${namePrefix}_${length}.log 2>/mnt/scratch/Homa/ChoppedReads/aligned_rCandy_${namePrefix}_${length}.err &
			#time ~/.cabal/bin/r-candy -o /mnt/scratch/Homa/ChoppedReads/aligned_rCandy_${namePrefix}_${length}.bam -A ${AS}  -G /r1/people/homa_amini/Repository/SimulateGenome/HighQualityReads.tmp.anfoIndex  /mnt/scratch/Homa/ChoppedReads/${namePrefix}_${length}.bam 1>/mnt/scratch/Homa/ChoppedReads/aligned_rCandy_${namePrefix}_${length}.log 2>/mnt/scratch/Homa/ChoppedReads/aligned_rCandy_${namePrefix}_${length}.err &
			#Mapped with r-candy without parameters
#			Old path /r1/people/homa_amini/.cabal/bin/r-candy
#			r-candy -o /mnt/scratch/Homa/ChoppedReads/aligned_rCandy_${namePrefix}_${length}_np.bam -A ${AS} -G /r1/people/homa_amini/Repository/SimulateGenome/SimulateGenome.tmp.anfoIndex /mnt/scratch/Homa/ChoppedReads/${namePrefix}_${length}.bam &&
			#Do some counting (e.g Match,Mismatch,w/o MAQ) for ploting
#			/home/homa_amini/Repository/Hash-rCandy/Hash_rCandy -i "/mnt/scratch/Homa/ChoppedReads/aligned_BWA_${namePrefix}_${length}.bam" -o "/mnt/scratch/Homa/ChoppedReads/BWA_" -m "BWA" -l ${length} 1 >> SimulateGenome.all.stat.tex
#			/home/homa_amini/Repository/Hash-rCandy/Hash_rCandy -i "/mnt/scratch/Homa/ChoppedReads/aligned_rCandy_${namePrefix}_${length}.bam" -o "/mnt/scratch/Homa/ChoppedReads/rCandy_" -m "rCandy" -l ${length} 1 >> SimulateGenome.all.stat.tex
#			/home/homa_amini/Repository/Hash-rCandy/Hash_rCandy -i "/mnt/scratch/Homa/ChoppedReads/aligned_rCandy_${namePrefix}_${length}_np.bam" -o "/mnt/scratch/Homa/ChoppedReads/rCandy_np_" -m "rCandy-np" -l ${length} 1 >> SimulateGenome.all.stat.tex
			# Just rate FPR/TPR
#			/home/homa_amini/Repository/Hash-rCandy/Hash_rCandy -i "/mnt/scratch/Homa/ChoppedReads/aligned_rCandy_${namePrefix}_${length}.bam" -o "/mnt/scratch/Homa/ChoppedReads/rCandy_${AS}_" -m "rCandy-${namePrefix}-${AS}" -l ${length} -r 1 >> SimulateGenome.all.Rate
#			/home/homa_amini/Repository/Hash-rCandy/Hash_rCandy -i "/mnt/scratch/Homa/ChoppedReads/aligned_rCandy_${namePrefix}_${length}_np.bam" -o "/mnt/scratch/Homa/ChoppedReads/rCandy_np_${AS}_" -m "rCandy-np-${namePrefix}-${AS}" -l ${length} -r 1 >> SimulateGenome.all.Rate
#			echo " \\hline" >> SimulateGenome.all.stat.tex


		done

	done
done
#echo "	\\end{longtable}" >> SimulateGenome.all.stat.tex
#echo "\\end{document}" >> SimulateGenome.all.stat.tex
#Draw Plot From R
#/usr/local64/bin/R --no-save --no-restore --args "/r1/people/homa_amini/Repository/SimulateGenome/fakeGenome.all.stat" "/r1/people/homa_amini/Repository/SimulateGenome/Rplot_SimulateGenome.pdf" < plot.R
#Make pdf from lates
#pdflatex SimulateGenome.all.stat.tex
#pdflatex SimulateGenome.all.stat.tex
#echo "Done"
