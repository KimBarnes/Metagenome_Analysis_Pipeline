directory=""
python $directory/Cluster_Analysis/FileSize_Filter.py
checkm lineage_wf -t 64 -x fasta $directory/ANALYSIS/ $directory/ANALYSIS/CheckM
python $directory/Cluster_Analysis/CheckM_Parser1.py
python $directory/Cluster_Analysis/Multi-Contig_Contamination_Removal.py
for D in $directory/ANALYSIS/Sequence_Composition/*;
do checkm lineage_wf -t 64 -x fasta $D/ $D/CheckM
done
python $directory/Cluster_Analysis/CheckM_Parser2.py
setup_prokka
cd ~/$directory/ANALYSIS/Genome_Annotation/Bacteria
for file in *.fasta;
do prokka -metagenome "$file" -outdir Prokka/$file -locustag $file -kingdom Bacteria -centre X -compliant
done
cd ~/$directory/ANALYSIS/Genome_Annotation/Archaea
for file in *.fasta;
do prokka -metagenome "$file" -outdir Prokka/$file -locustag $file -kingdom Archaea -centre X -compliant
done
cd ~
python $directory/Cluster_Analysis/Metabolic_Pathway_Prediction.py
