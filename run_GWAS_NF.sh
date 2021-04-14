## Load modules
ml purge
ml nextflow/21.03.0-edge
ml Singularity/3.4.2

WORDIR=/camp/stp/babs/scratch/schneid/Regenie_GWA_NF/

mkdir -p ${DIR}

nextflow Regenie_GWAS_NF.nf --input -work-dir $WORKDIR -resume
