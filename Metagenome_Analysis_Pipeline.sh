python Final_Test/Cluster_Analysis/FileSize.py
checkm lineage_wf -t 64 -x fasta Final_Test/ANALYSIS/ Final_Test/ANALYSIS/CheckM
python Final_Test/Cluster_Analysis/CheckM_Reader.py
python Final_Test/Cluster_Analysis/Genome_Annotation.py
python Final_Test/Cluster_Analysis/Sequence_Analysis.py
for D in Final_Test/ANALYSIS/Sequence_Composition/*;
do checkm lineage_wf -t 64 -x fasta $D/ $D/CheckM
done
python Final_Test/Cluster_Analysis/Genome_Annotation2.py
setup_prokka
cd ~/Final_Test/ANALYSIS/Genome_Annotation/Bacteria
for file in *.fasta;
do prokka -metagenome "$file" -outdir Prokka/$file -locustag $file -kingdom Bacteria -centre X -compliant
done
cd ~/Final_Test/ANALYSIS/Genome_Annotation/Archaea
for file in *.fasta;
do prokka -metagenome "$file" -outdir Prokka/$file -locustag $file -kingdom Archaea -centre X -compliant
done
cd ~
python Final_Test/Cluster_Analysis/Metabolic_Pathways_Linux.py
