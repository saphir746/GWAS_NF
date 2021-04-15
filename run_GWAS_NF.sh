## Load modules
ml purge
ml nextflow/21.03.0-edge
ml Singularity/3.4.2

WORDIR=/camp/stp/babs/scratch/schneid/Regenie_GWA_NF/

mkdir -p ${WORKDIR}

data_dir="/camp/stp/babs/working/schneid/projects/swantonc/aska.przewrocka/ZNF516_UKB/Data_DSL/UKB_ZNF516/RegenieXNextFlow/Regenie_PhenoPre/results/"
genom_dir="/camp/stp/babs/working/schneid/projects/swantonc/aska.przewrocka/ZNF516_UKB/Data_DSL/UKB_ZNF516/RegenieXNextFlow/Regenie_NF_main/results_bed/"

nextflow Regenie_GWAS_NF.nf -resume \
	--Cat_data ${data_dir}"ZNF516_CHD.txt" \
	--Cont_data ${data_dir}"ZNF516_quant_traits.txt" \
	--Covar_data ${data_dir}"ZNF516_covariates.txt" \
	--Excl_ids ${data_dir}"ZNF516_genetic_excl.csv" \
	--sample "/camp/project/proj-tracerx-lung/public_datasets/UK_BIOBANK/ImputedGenomicData/ukb51637_imp_chr18_v3_s487298.sample" \
	--out_ch "/camp/stp/babs/working/schneid/projects/swantonc/aska.przewrocka/ZNF516_UKB/Data_DSL/UKB_ZNF516/RegenieXNextFlow/Regenie_NF_main/results/" \
        --Bed ${genom_dir}"ukb_cal_Chr18.bed" \
	--Bim ${genom_dir}"ukb_cal_Chr18.bim" \
	--Fam ${genom_dir}"ukb_cal_Chr18.fam" \
	--list1 ${genom_dir}"ukb_cal_step1.snplist" \
	-work-dir $WORKDIR 
