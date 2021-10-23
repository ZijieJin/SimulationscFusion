# SimulationscFusion

This project generates simulated single-cell RNA-seq files with given gene fusions.

# Generating Reads from Given Fusion Transcripts

We use [FusionSimulatorToolkit](https://github.com/FusionSimulatorToolkit/FusionSimulatorToolkit/wiki) to generate reads from fusion transcripts.

The following software are required in addition to the scripts provided here:

* [RSEM](http://deweylab.github.io/RSEM/)
* [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
* [Bowtie-1](http://bowtie-bio.sourceforge.net/index.shtml)

Run scripts below:

`$TRINITY_HOME/util/align_and_estimate_abundance.pl 
                    --transcripts TrueFusion.fa 
                    --est_method RSEM 
                    --aln_method bowtie 
                    --prep_reference 
                    --seqType fq 
                    --left adipose-ERR030880_1.fastq.gz.revised.gz 
                    --right adipose-ERR030880_2.fastq.gz.revised.gz 
                    --output_dir RSEM_outdir`

`simulate_fusion_trans_expr_vals.pl 
           RSEM_outdir/RSEM.isoforms.results 
           TrueFusion.fa  
        > target.forSimulation.RSEM.isoforms.results`

`$RSEM_HOME/rsem-simulate-reads TrueFusion.fa.RSEM
                              RSEM_outdir/RSEM.stat/RSEM.model
                              target.forSimulation.RSEM.isoforms.result
                              0.01  2000  sim_reads`


# Generating Background Reads

`$TRINITY_HOME/util/align_and_estimate_abundance.pl 
                    --transcripts hg19.cdna.fa 
                    --est_method RSEM 
                    --aln_method bowtie 
                    --prep_reference 
                    --seqType fq 
                    --left adipose-ERR030880_1.fastq.gz.revised.gz 
                    --right adipose-ERR030880_2.fastq.gz.revised.gz 
                    --output_dir RSEM_outdir`

`simulate_fusion_trans_expr_vals.pl 
           RSEM_outdir/RSEM.isoforms.results 
           hg19.cdna.fa  
        > target.forSimulation.RSEM.isoforms.results`

`$RSEM_HOME/rsem-simulate-reads hg19.cdna.fa.RSEM
                              RSEM_outdir/RSEM.stat/RSEM.model
                              target.forSimulation.RSEM.isoforms.result
                              0.01  4000000  sim_reads`


