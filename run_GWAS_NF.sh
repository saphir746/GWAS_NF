## Load modules
ml purge
ml Nextflow/20.12.0-edge
ml Graphviz/2.38.0-foss-2016b
## mlSingularity/3.4.2

WORDIR=/camp/stp/babs/scratch/schneid/Regenie_GWA_NF/
OUTDIR=/camp/stp/babs/working/schneid/projects/swantonc/aska.przewrocka/ZNF516_UKB/Data_DSL/UKB_ZNF516/GWAS_NF/results

if [ ! -d $WORKDIR ]
then
	mkdir $WORKDIR
fi

if [ ! -d $OUTDIR ]
then
	mkdir $OUTDIR
fi

UKB_dir="/camp/project/proj-tracerx-lung/public_datasets/UK_BIOBANK/ImputedGenomicData/"
data_dir="/camp/stp/babs/working/schneid/projects/swantonc/aska.przewrocka/ZNF516_UKB/Data_DSL/UKB_ZNF516/RegenieXNextFlow/Regenie_PhenoPre/results/"
genom_dir="/camp/stp/babs/working/schneid/projects/swantonc/aska.przewrocka/ZNF516_UKB/Data_DSL/UKB_ZNF516/RegenieXNextFlow/Regenie_NF_main/results_bed/"

nextflow Regenie_GWAS_NF.nf -with-dag GWAS.svg -resume \
	--Cat_data ${data_dir}"ZNF516_CHD.txt" \
	--Cont_data ${data_dir}"ZNF516_quant_traits.txt" \
	--Covar_data ${data_dir}"ZNF516_covariates.txt" \
	--Excl_ids ${data_dir}"ZNF516_genetic_excl.csv" \
	--sample ${UKB_dir}"ukb51637_imp_chr18_v3_s487298.sample" \
	--bgen ${UKB_dir}"ukb_imp_chr18_v3.bgen" \
	--out_ch $OUTDIR \
	-work-dir "/camp/stp/babs/scratch/schneid/Regenie_GWA_NF/" \
        --Bed ${genom_dir}"ukb_cal_Chr18.bed" \
	--Bim ${genom_dir}"ukb_cal_Chr18.bim" \
	--Fam ${genom_dir}"ukb_cal_Chr18.fam" \
	--list1 ${genom_dir}"ukb_cal_step1.snplist"
