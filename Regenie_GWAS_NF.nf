#!/usr/bin/env nextflow

import java.nio.file.Paths

// scripts
project_dir=Paths.get(System.getProperty("user.dir"),"scripts/").toString()
excl_script = Paths.get(project_dir, "genom_exclusion.py")
manPlot_script = Paths.get(project_dir, "GWA_results.py")
GWA_post = Paths.get(project_dir, "GWA_postprocess.py")

/*
* * 1.B: Sample exclusions
*/
process genom_excl {

    executor = "local"

    conda '/camp/stp/babs/working/software/anaconda/envs/biobankread/'

    output:
    file 'qc_pass2.id' into Excl_ch_2

    script:
    """
    #!/usr/bin/env bash
    
    python ${excl_script} --genomExcl ${params.Excl_ids}
    """
}

Excl_ch_2.into{ Excl_ch_bt1; Excl_ch_bt2; Excl_ch_cont1; Excl_ch_cont2 }

/*
* * 2: Regenie - step 1
* *
*/
  
process Regenie_1 { 
 
    conda '/camp/stp/babs/working/schneid/conda/envs/RegenieGWA'
 
    input:
    file exlcIDs from Excl_ch_bt1
 
    output:
    file('ukb_step1_BT_pred.list') into regenie1_ch_bt
 
    script:
    """
    ln -s ${params.Bed} ${params.Bim} ${params.Fam} .

    regenie --step 1 \
        --bed ukb_cal_Chr18 \
        --extract ${params.list1} \
        --keep $exlcIDs \
        --covarFile ${params.Covar_data} \
        --covarColList "Age","Gender" \
        --phenoFile ${params.Cat_data} \
        --bt \
        --bsize 701 \
        --out ukb_step1_BT
    """
}


process Regenie_1_cont {

    conda '/camp/stp/babs/working/schneid/conda/envs/RegenieGWA'

    input:
    file exlcIDs from Excl_ch_cont1

    output:
    file('ukb_step1_cont_pred.list') into regenie1_ch_cont

    script:
    """
    ln -s ${params.Bed} ${params.Bim} ${params.Fam} .

    regenie --step 1 \
        --bed ukb_cal_Chr18 \
        --extract ${params.list1} \
        --keep $exlcIDs \
        --covarFile ${params.Covar_data} \
        --covarColList "Age","Gender" \
        --phenoFile ${params.Cont_data} \
        --bsize 701 \
        --lowmem \
        --out ukb_step1_cont
    """
}


/**
*** 3: Regenie - step 2
***/

process Regenie_2 {

    errorStrategy 'ignore'
    conda '/camp/stp/babs/working/schneid/conda/envs/RegenieGWA'

    input:
    file('ukb_step1_BT_pred.list') from regenie1_ch_bt
    file exlcIDs from Excl_ch_bt2

    output:
    file('*.regenie') into regenie_ch_final

    script:
    """
    regenie --step 2 \
        --bgen ${params.bgen} \
	--ref-first \
	--sample ${params.sample} \
        --keep $exlcIDs \
        --pred ukb_step1_BT_pred.list \
        --phenoFile ${params.Cat_data} \
        --covarFile ${params.Covar_data} \
        --covarColList "Age","Gender" \
        --bt \
        --firth --approx \
        --bsize 500 \
        --split \
        --out ukb_step2_BT_chr18
    """
}

//process Regenie_2_cont {
//
//    errorStrategy 'ignore'
//    conda '/camp/stp/babs/working/schneid/conda/envs/RegenieGWA'
//
//    input:
//    file('ukb_step1_cont_pred.list') from regenie1_ch_cont
//    file exlcIDs from Excl_ch_cont2
//
//    output:
//    file('*.regenie') into regenie_ch_final_cont
//
//    script:
//    """
//    ln -s ${params.Bed} ${params.Bim} ${params.Fam} .
//
//    regenie --step 2 \
//        --bed ukb_cal_Chr18 \
//        --extract ${params.list1} \
//        --keep $exlcIDs \
//        --pred ukb_step1_cont_pred.list \
//        --phenoFile ${params.Cont_data} \
//        --covarFile ${params.Covar_data} \
//        --covarColList "Age","Gender" \
//        --bsize 500 \
//        --spa --approx \
//        --split \
//        --out ukb_step2_chr18
//    """
//}


process Regenie_2_cont {

    errorStrategy 'ignore'
    label 'process_high'
    conda '/camp/stp/babs/working/schneid/conda/envs/RegenieGWA'

    input:
    file('ukb_step1_cont_pred.list') from regenie1_ch_cont
    file exlcIDs from Excl_ch_cont2

    output:
    file('*.regenie') into regenie_ch_final_cont

    script:
    """
    regenie --step 2 \
        --bgen ${params.bgen} \
        --ref-first \
	--sample ${params.sample} \
        --pred ukb_step1_cont_pred.list \
	--keep $exlcIDs \
        --phenoFile ${params.Cont_data} \
        --covarFile ${params.Covar_data} \
        --covarColList "Age","Gender" \
        --bsize 500 \
        --spa --approx \
        --split \
        --out ukb_step2_chr18
    """
}


regenie_ch_final.flatten().concat(regenie_ch_final_cont.flatten()).set{ regenie_ch_final_all }

process GWA_summary {

   executor 'local'


   input:
   file regenie from regenie_ch_final_all

   output:
   file('ZNF516_GWA_*_results.csv') into Ch_final
   val 'OK' into Ch_trigger

   publishDir "${params.out_ch}/", mode: params.publish_dir_mode, overwrite: true,
   saveAs: { filename -> if ( filename.indexOf(".csv") != 1 ) {
                                     "$filename"
                              } else {
                                      null
                              }
      }

   script:
   """
   #!/usr/bin/env bash
   
   python ${manPlot_script} --RegenieRes ${regenie}
   """
}                  

process GWA_postProcess {
   executor 'local'

   conda '/camp/stp/babs/working/software/anaconda/envs/biobankread'

   input:
   val 'OK' from Ch_trigger

   output:
   file('ZNF516_All_GWA_signif_results.csv') into Ch_signif

   publishDir "${params.out_ch}/", mode: params.publish_dir_mode, overwrite: true,
   saveAs: { filename -> if ( filename.indexOf(".csv") != 1 ) {
                                        "$filename"
                                } else {
                                        null
                                }
        }

   script:
   """
   #!/usr/bin/env bash
   
   python ${GWA_post} --Folder ${params.out_ch}/
   """
}

workflow.onComplete {
}         
