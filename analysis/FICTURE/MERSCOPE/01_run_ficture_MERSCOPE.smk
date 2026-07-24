## run ficture scripts for transcripts in ROIs
## Per sample based

sample_names =[
    	"23MH0007",
    	"23MH0026"
    ]

transcript_dir = "transcripts_qc/"
outdir = "predefine_model_m/"
ficture_dir = "../"
def get_mb(wildcards,attempt):
    # This function can be implemented to dynamically determine the memory requirement based on the sample or other criteria.
    # For now, it returns a hardcoded value of 120000 MB (120 GB).
    return 120000 + attempt*50000
rule all:
    input:
#        expand(outdir + "data/{s_n}/analysis/nF20.d_12/figure/nF20.d_12.decode.prj_12.r_4_5.pixel.png", s_n=sample_names),
        # expand(outdir + "data/{s_n}/analysis/{model_id}/figure/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel.png", s_n=sample_names, model_id='nF20.d_12', fit_width='12', anchor_res='4', radius='5'),
        # expand(outdir + "data/{s_n}/analysis/{model_id}/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel.tsv.gz", s_n=sample_names, model_id='nF20.d_12', fit_width='12', anchor_res='4', radius='5'),
        # expand( outdir + "data/{s_n}/analysis/{model_id}/figure/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel_set1.png", s_n='KLH_417', model_id='nF20.d_12', fit_width='12', anchor_res='4', radius='5'),   
        expand( outdir + "data/{s_n}/analysis/{model_id}/figure/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel_{set}.png", s_n=sample_names, model_id='nF20.d_12', fit_width='12', anchor_res='4', radius='5',
                                                                                                                                            set=['All']),
        expand( outdir + "data/{s_n}/analysis/{model_id}/figure/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel_noncat_{set}.png", s_n=sample_names, model_id='nF20.d_12', fit_width='12', anchor_res='4', radius='5',
                                                                                                                                            set=['All'])     

rule prep_input_tsv:
    input:
        roi_tx = transcript_dir + "{s_n}/{s_n}.csv",
        script = ficture_dir + "ficture_scripts/00_prep_input.sh"
    output:
        tsv = outdir + "data/{s_n}/{s_n}.tsv.gz"
    resources:
        mem=16000,
        time="02:00:00",
        cpus=1
    params:
        slurm_err = outdir + "data/{s_n}/prep_input.err",
        slurm_out = outdir + "data/{s_n}/prep_input.out",
    shell:
        """
            sh {input.script} {input.roi_tx} {output.tsv}
        """

rule get_min_max:
    input:
        tsv = outdir + "data/{s_n}/{s_n}.tsv.gz",
        script = ficture_dir + "ficture_scripts/01_get_min_max.sh"
    output:
        minmax = outdir + "data/{s_n}/coordinate_minmax.tsv"
    resources:
        mem=16000,
        time="02:00:00",
        cpus=1
    params:
        slurm_err = outdir + "data/{s_n}/get_min_max.err",
        slurm_out = outdir + "data/{s_n}/get_min_max.out"
    shell:
        """
            sh {input.script} {input.tsv} {output.minmax}
        """

rule run_prep_minibatch:
    input: 
        tsv = outdir + "data/{s_n}/{s_n}.tsv.gz",
        coord_tsv = outdir + "data/{s_n}/coordinate_minmax.tsv",
        script = ficture_dir + "ficture_scripts/generic_submit.sh"
    params:
        path=outdir + "data/{s_n}",
        slurm_out = outdir + "data/{s_n}/ficture_prep_minibatch.out",
        slurm_err = outdir + "data/{s_n}/ficture_prep_minibatch.err",
        mu_scale = 1,
        major_axis ='X',
        batch_size=1000,
        batch_buff=60
    resources:
        mem=86000,
        time="2:00:00",
        cpus=4
    output:
        batch = outdir + "data/{s_n}/batched.matrix.tsv.gz"
        
    shell:
        """
            module load micromamba/ 
            eval "$(micromamba shell hook --shell bash)"
            micromamba activate /vast/projects/CDFC/rlyu/Softwares/condaenvs/ficture
            echo {wildcards.s_n}
            # Prepare training minibatches, only need to run once if you plan to fit multiple models (say with different number of factors)
            batch_gz={output.batch}
            batch=$(echo $batch_gz | sed 's/\.gz$//g')
 
            ficture make_spatial_minibatch \
                        --input {input.tsv} \
                        --output $batch \
                        --mu_scale {params.mu_scale} \
                        --batch_size {params.batch_size} \
                        --batch_buff {params.batch_buff} \
                        --major_axis {params.major_axis}
            
            
            gzip $batch


        """
rule run_transform:
    input: 
        batch = outdir + "data/{s_n}/batched.matrix.tsv.gz",
        tsv = outdir + "data/{s_n}/{s_n}.tsv.gz",

        coord_tsv = outdir + "data/{s_n}/coordinate_minmax.tsv",
        model = outdir + "model_matrix/{s_n}/model_matrix.tsv.gz",
        script = ficture_dir + "ficture_scripts/generic_submit.sh"
    params:
        path=outdir + "data/{s_n}",
        slurm_out = outdir + "data/{s_n}/ficture_prep_minibatch.out",
        slurm_err = outdir + "data/{s_n}/ficture_prep_minibatch.err",
        mu_scale = 1,
        major_axis ='X',
        batch_size=1000,
        batch_buff=60,
        key='Count',
        min_ct_per_unit=50,
        min_ct_per_unit_fit=20,
        min_ct_per_feature=50,
        R=10,
        fit_width= '12',
        fit_nmove='3',#((12 / anchor_res))
        anchor_res=4,
        seed=20251,
        model_id='nF20.d_12',
        # For DE output
        max_pval_output=1e-3,
        min_fold_output=1.5,
    resources:
        mem=86000,
        time="2:00:00",
        cpus=4
    output:
        anchor_fit_results = outdir + "data/{s_n}/analysis/{model_id}/{model_id}.prj_{fit_width}.r_{anchor_res}.fit_result.tsv.gz"
    shell:
        """
            module load micromamba/ 
            eval "$(micromamba shell hook --shell bash)"
            micromamba activate /vast/projects/CDFC/rlyu/Softwares/condaenvs/ficture
            # Coarse plot for inspection
            
            output_path={params.path}/analysis/{wildcards.model_id}
            anchor_info=prj_{wildcards.fit_width}.r_{wildcards.anchor_res}

            # Transform
            output=$output_path/{wildcards.model_id}.$anchor_info
            ficture transform --input {input.tsv} \
                             --output_pref $output --model {input.model} \
                             --key {params.key} --major_axis {params.major_axis} --hex_width {wildcards.fit_width} \
                             --n_move {params.fit_nmove} --min_ct_per_unit {params.min_ct_per_unit_fit} \
                             --mu_scale {params.mu_scale} --thread {resources.cpus} --precision 2
        """
# + ficture transform --input ../../output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/data/WT_37RF/WT_37RF.tsv.gz
#  --output_pref ../../output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/data/WT_37RF/
#  --model /vast/projects/CDFC/rlyu/Projects/Mouse_Myeloma_2025/output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/model_matrix.tsv.gz 
#  --key Count --major_axis X --hex_width 12 --n_move 3 --min_ct_per_unit 20 --mu_scale 1 --thread 16 --precision 2


# rule run_ficture_decode:
#     input: 
#         tsv = outdir + "data/{s_n}/{s_n}.tsv.gz",
#         coord_tsv = outdir + "data/{s_n}/coordinate_minmax.tsv",
#         script = ficture_dir + "ficture_scripts/05_run_prep_minibatch_and_decode.sh"
#     params:
#         path=outdir + "data/{s_n}"
#     resources:
#         mem=86000,
#         time="12:00:00",
#         cpus=16
#     output:
#         ficture_out = outdir + "data/{s_n}/analysis/nF20.d_12/figure/nF20.d_12.decode.prj_12.r_4_5.pixel.png"
#     shell:
#         """
#             sh {input.script} {params.path} {wildcards.s_n} 
#         """

# radius = $(($anchor_res+1))
rule run_ficture_decode:
    input: 
        pixel = outdir + "data/{s_n}/batched.matrix.tsv.gz",
        model = outdir + "model_matrix/{s_n}/model_matrix.tsv.gz",
        anchor_fit_results = outdir + "data/{s_n}/analysis/{model_id}/{model_id}.prj_{fit_width}.r_{anchor_res}.fit_result.tsv.gz",
        cmap= outdir + "model_matrix/{s_n}/celltype.rgb.tsv"
    params:
        slurm_err = outdir + "data/{s_n}/ficture_decode.err",
        slurm_out = outdir + "data/{s_n}/ficture_decode.out",
        radius=5,
        mu_scale=1,
        key='Count',
        precision=0.1,
        topk=3,
        output_prefix=outdir + "data/{s_n}/analysis/{model_id}/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}",
        min_ct_per_feature='50', 
        max_pval_output='1e-3',
        min_fold_output='1.5',
        path= outdir + "data/{s_n}/analysis/{model_id}",
        bulk_chisq = outdir + "data/{s_n}/analysis/{model_id}/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.bulk_chisq.tsv"
    resources:
        mem=get_mb,
        time="12:00:00",
        cpus=16
    output:
        decode_out = outdir + "data/{s_n}/analysis/{model_id}/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel.tsv.gz",
        posterior_count = outdir + "data/{s_n}/analysis/{model_id}/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.posterior.count.tsv.gz"
    shell:
        """
            module load micromamba/ 
            eval "$(micromamba shell hook --shell bash)"
            micromamba activate /vast/projects/CDFC/rlyu/Softwares/condaenvs/ficture       

            ficture slda_decode --input {input.pixel} --output {params.output_prefix} \
                                --model {input.model} --anchor {input.anchor_fit_results} \
                                --anchor_in_um --neighbor_radius {wildcards.radius} --mu_scale {params.mu_scale} \
                                --key {params.key} --precision {params.precision} --lite_topk_output_pixel {params.topk} \
                                --lite_topk_output_anchor {params.topk} --thread {resources.cpus}
                                
            # Ficture factor feature report
            ficture de_bulk --input {output.posterior_count} --output {params.bulk_chisq} --min_ct_per_feature {params.min_ct_per_feature} \
                            --max_pval_output {params.max_pval_output} --min_fold_output {params.min_fold_output} \
                            --thread {resources.cpus}
            ficture factor_report --path {params.path} \
                                  --pref {wildcards.model_id}.decode.prj_{wildcards.fit_width}.r_{wildcards.anchor_res}_{wildcards.radius} \
                                  --color_table {input.cmap}

        """

# + ficture slda_decode 
# --input ../../output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/data/WT_37RF/batched.matrix.tsv.gz 
# --output ../../output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/data/WT_37RF/analysis/nF20.d_12/nF20.d_12.decode.prj_12.r_4_5 
# --model /vast/projects/CDFC/rlyu/Projects/Mouse_Myeloma_2025/output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/model_matrix.tsv.gz 
# --anchor ../../output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/data/WT_37RF/analysis/nF20.d_12/nF20.d_12.prj_12.r_4.fit_result.tsv.gz 
# --anchor_in_um 
# --neighbor_radius 5 --mu_scale 1 --key Count --precision 0.1 --lite_topk_output_pixel 3 --lite_topk_output_anchor 3 --thread 16

rule sort_and_index_pixel_result:
    input:
        pixel_tsv = outdir + "data/{s_n}/analysis/{model_id}/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel.tsv.gz",
        coordinate_minmax = outdir + "data/{s_n}/coordinate_minmax.tsv"
    output:
        sorted_pixel_tsv =outdir + "data/{s_n}/analysis/{model_id}/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel.sorted.tsv.gz"
    params:
        slurm_err = outdir + "data/{s_n}/ficture_sort_pixel.err",
        slurm_out = outdir + "data/{s_n}/ficture_sort_pixel.out",
        bsize = 2000,
        scale = 100,  # Replace with actual value or pass dynamically
        K = 3  # Replace with actual value or pass dynamically
    resources:
        mem = 30000,
        cpus=1
    shell:
        """
            module load micromamba/ 
            eval "$(micromamba shell hook --shell bash)"
            micromamba activate /vast/projects/CDFC/rlyu/Softwares/condaenvs/ficture  
            # Extract xmin, ymin, xmax, ymax from the coordinate_minmax file
            xmin=$(awk '$1 == "xmin" {{print $2}}' {input.coordinate_minmax})
            xmax=$(awk '$1 == "xmax" {{print $2}}' {input.coordinate_minmax})
            ymin=$(awk '$1 == "ymin" {{print $2}}' {input.coordinate_minmax})
            ymax=$(awk '$1 == "ymax" {{print $2}}' {input.coordinate_minmax})

            # Define variables
            input={input.pixel_tsv}
            output={output.sorted_pixel_tsv}
            offsetx=${{xmin}}
            offsety=${{ymin}}
            bsize={params.bsize}
            scale={params.scale}
            rangex=$( echo "(${{xmax}} - ${{xmin}} + 0.5)/1+1" | bc )
            rangey=$( echo "(${{ymax}} - ${{ymin}} + 0.5)/1+1" | bc )

            # Define header
            header="##K={params.K};TOPK=3\\n##BLOCK_SIZE={params.bsize};BLOCK_AXIS=X;INDEX_AXIS=Y\\n##OFFSET_X=${{offsetx}};OFFSET_Y=${{offsety}};SIZE_X=${{rangex}};SIZE_Y=${{rangey}};SCALE={params.scale}\\n#BLOCK\\tX\\tY\\tK1\\tK2\\tK3\\tP1\\tP2\\tP3"

            # Process and sort the file
            (echo -e "${{header}}" && zcat ${{input}} | tail -n +2 | perl -slane '$F[0]=int(($F[1]-$offx)/$bsize) * $bsize; $F[1]=int(($F[1]-$offx)*$scale); $F[1]=($F[1]>=0)?$F[1]:0; $F[2]=int(($F[2]-$offy)*$scale); $F[2]=($F[2]>=0)?$F[2]:0; print join("\\t", @F);' -- -bsize={params.bsize} -scale={params.scale} -offx=${{offsetx}} -offy=${{offsety}} | sort -S {resources.mem}M -k1,1g -k3,3g ) | gzip -c > ${{output}}

        """

rule de_report:
    input:
        posterior_count = outdir + "data/{s_n}/analysis/{model_id}/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.posterior.count.tsv.gz",
        cmap= outdir + "model_matrix/{s_n}/celltype.rgb.tsv"
    params:
        min_ct_per_feature='50', 
        max_pval_output='1e-3' ,
        min_fold_output='1.5',
        path= outdir + "data/{s_n}/analysis/{model_id}",
        slurm_err= outdir + "data/{s_n}/ficture_de_report.err",
        slurm_out=outdir + "data/{s_n}/ficture_de_report.out",
    resources:
        cpus=16,
        mem=20000

    output:
        bulk_chisq = outdir + "data/{s_n}/analysis/{model_id}/DE/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.bulk_chisq.tsv",
    shell:
        """
            module load micromamba/ 
            eval "$(micromamba shell hook --shell bash)"
            micromamba activate /vast/projects/CDFC/rlyu/Softwares/condaenvs/ficture  

            ficture de_bulk --input {input.posterior_count} --output {output.bulk_chisq} --min_ct_per_feature ${min_ct_per_feature} \
                            --max_pval_output {params.max_pval_output} --min_fold_output {params.min_fold_output} \
                            --thread {resources.cpus}

            output=${output_path}/${prefix}.factor.info.html
            ficture factor_report --path {params.path} \
                                  --pref {widlcards.model_id}.decode.prj_{widlcards.fit_width}.r_{widlcards.anchor_res}_{widlcards.radius} \
                                  --color_table {input.cmap}
        """
# ficture de_bulk --input ${input} --output ${output} --min_ct_per_feature ${min_ct_per_feature} --max_pval_output ${max_pval_output} --min_fold_output ${min_fold_output} --thread ${thread}
# + ficture de_bulk 
# --input ../../output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/data/WT_37RF/analysis/nF20.d_12/nF20.d_12.decode.prj_12.r_4_5.posterior.count.tsv.gz 
# --output ../../output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/data/WT_37RF/analysis/nF20.d_12/DE/nF20.d_12.decode.prj_12.r_4_5.bulk_chisq.tsv 
# --min_ct_per_feature 50 --max_pval_output 1e-3 --min_fold_output 1.5 --thread 16

# + ficture factor_report 
# --path ../../output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/data/WT_37RF/analysis/nF20.d_12 
# --pref nF20.d_12.decode.prj_12.r_4_5 \
# --color_table /vast/projects/CDFC/rlyu/Projects/Mouse_Myeloma_2025/output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/celltype.rgb.tsv
rule plot_full:
    input:
        sorted_pixel_tsv =outdir + "data/{s_n}/analysis/{model_id}/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel.sorted.tsv.gz",
        cmap= outdir + "model_matrix/{s_n}/celltype.rgb.tsv"
    params:
        pixel_resolution='0.5',
        path= outdir + "data/{s_n}/analysis/{model_id}",
        slurm_err= outdir + "data/{s_n}/ficture_plot_pixel.err",
        slurm_out=outdir + "data/{s_n}/ficture_plot_pixel.out"
    resources:
        cpus=16,
        mem=20000

    output:
        png_full = outdir + "data/{s_n}/analysis/{model_id}/figure/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel.png",
    shell:
        """
            module load micromamba/ 
            eval "$(micromamba shell hook --shell bash)"
            micromamba activate /vast/projects/CDFC/rlyu/Softwares/condaenvs/ficture  

            ficture plot_pixel_full --input {input.sorted_pixel_tsv} --color_table {input.cmap} \
                                    --output {output.png_full} --plot_um_per_pixel {params.pixel_resolution} --full

        """
# + ficture plot_pixel_full 
# \--input ../../output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/data/WT_37RF/analysis/nF20.d_12/nF20.d_12.decode.prj_12.r_4_5.pixel.sorted.tsv
# .gz
#  --color_table /vast/projects/CDFC/rlyu/Projects/Mouse_Myeloma_2025/output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/celltype.rgb.tsv 
#  --output ../../output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/data/WT_37RF/analysis/nF20.d_12/figure/nF20.d_12.decode.prj_12.r_4_5.pixel.png --plot_um_per_pixel 0.5 --full


rule plot_pixel:
    input:
        sorted_pixel_tsv =outdir + "data/{s_n}/analysis/{model_id}/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel.sorted.tsv.gz",
        cmap= outdir + "model_matrix/{s_n}/celltype.rgb.tsv",
        script = ficture_dir + "ficture_scripts/plot_pixel.py"
    params:
        pixel_resolution='0.5',
        path= outdir + "data/{s_n}/analysis/{model_id}",
        slurm_err= outdir + "data/{s_n}/plot_pixel_highlight.err",
        slurm_out=outdir + "data/{s_n}/plot_pixel_highlight.out"
    resources:
        cpus=2,
        mem=120000

    output:
        png_highlight = outdir + "data/{s_n}/analysis/{model_id}/figure/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel_set1.png"
    shell:
        """
            module load micromamba/ 
            eval "$(micromamba shell hook --shell bash)"
            micromamba activate /vast/projects/CDFC/rlyu/Softwares/condaenvs/spatialrna  

            python {input.script} \
            --input_tsv  {input.sorted_pixel_tsv} \
            --cmap_dict_tsv {input.cmap} \
            --output {output.png_highlight} \
            --pixel_size 50 # \
#             --highlight_clusters 0 1 21 11 12 3 8  

        """
def get_channels_to_highlight(wildcards):
    # This function can be implemented to dynamically determine which channels to highlight based on the data or other criteria.
    # For now, it returns a hardcoded list of channels.
    if wildcards.set == "All":
        return "-1"
    if wildcards.set == "Plasma":
        return "6"
    if wildcards.set == "Myeloid":
        return "5" # Myeloid
    if wildcards.set == "Stromal":
        return "8"  # Stromal
    if wildcards.set == "Immune": # B cell & Macrophage & Megs & Monocyte & T cell
        return "0 2 3 4 9"
    if wildcards.set == "HSC":
        return "1 7" # HSC & red cells
    return "0"
    
rule plot_pixel_multi:
    input:
        sorted_pixel_tsv =outdir + "data/{s_n}/analysis/{model_id}/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel.sorted.tsv.gz",
        cmap= outdir + "model_matrix/{s_n}/celltype.rgb.tsv",
    params:
        pixel_resolution='0.5',
        slurm_err= outdir + "data/{s_n}/plot_pixel_multi_{set}.err",
        slurm_out=outdir + "data/{s_n}/plot_pixel_multi_{set}.out",
        channel_list=get_channels_to_highlight
    resources:
        cpus=2,
        mem=120000
    wildcard_constraints:
        set="(?!noncat).*" 
    output:
        png_highlight = outdir + "data/{s_n}/analysis/{model_id}/figure/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel_{set}.png"
    shell:
        """
            module load micromamba/ 
            eval "$(micromamba shell hook --shell bash)"
            micromamba activate /vast/projects/CDFC/rlyu/Softwares/condaenvs/ficture  
            channels='{params.channel_list}'

            tmp_cmap=$(mktemp)

            awk -F'\t' -v chans="$channels" '
            BEGIN {{
                OFS="\t"
                
                if (chans == "-1") {{
                  keep_all = 1
                }} else {{
                  split(chans, arr, " ")
                  for (i in arr) keep[arr[i]] = 1
                }}
            }}
            NR==1 {{print; next}}

            # Keep highlighted factors
            (keep_all || ($1 in keep)) {{print; next}}

            # Fade everything else
            {{
                $3=0.45; $4=0.45; $5=0.45;
                print
            }}' {input.cmap} > $tmp_cmap
            cat $tmp_cmap 
            n=$(($(wc -l < "$tmp_cmap") - 1))
            channel_list=$(seq 0 $((n-1)) | tr '\n' ' ')
            echo $channel_list
            ficture plot_pixel_multi  --input {input.sorted_pixel_tsv} \
            --output {output.png_highlight} \
            --color_table $tmp_cmap \
            --plot_um_per_pixel {params.pixel_resolution} \
            --channel_list $channel_list   \
            --pcut 0.05 --categorical \
            --spcut 0.2

        """

rule plot_pixel_multi_noncategorical:
    input:
        sorted_pixel_tsv =outdir + "data/{s_n}/analysis/{model_id}/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel.sorted.tsv.gz",
        cmap= outdir + "model_matrix/{s_n}/celltype.rgb.tsv",
        script = ficture_dir + "ficture_scripts/plot_pixel.py"
    params:
        pixel_resolution='0.5',
        slurm_err= outdir + "data/{s_n}/plot_pixel_multi_{set}.err",
        slurm_out=outdir + "data/{s_n}/plot_pixel_multi_{set}.out",
        channel_list=get_channels_to_highlight
    resources:
        cpus=2,
        mem=120000

    output:
        png_highlight = outdir + "data/{s_n}/analysis/{model_id}/figure/{model_id}.decode.prj_{fit_width}.r_{anchor_res}_{radius}.pixel_noncat_{set}.png"
    shell:
        """
            module load micromamba/ 
            eval "$(micromamba shell hook --shell bash)"
            micromamba activate /vast/projects/CDFC/rlyu/Softwares/condaenvs/ficture  
            channels='{params.channel_list}'

            tmp_cmap=$(mktemp)

            awk -F'\t' -v chans="$channels" '
            BEGIN {{
                OFS="\t"
                
                if (chans == "-1") {{
                  keep_all = 1
                }} else {{
                  split(chans, arr, " ")
                  for (i in arr) keep[arr[i]] = 1
                }}
            }}
            NR==1 {{print; next}}

            # Keep highlighted factors
            (keep_all || ($1 in keep)) {{print; next}}

            # Fade everything else
            {{
                $3=0.45; $4=0.45; $5=0.45;
                print
            }}' {input.cmap} > $tmp_cmap
            cat $tmp_cmap 
            n=$(($(wc -l < "$tmp_cmap") - 1))
            channel_list=$(seq 0 $((n-1)) | tr '\n' ' ')
            echo $channel_list
            ficture plot_pixel_multi  --input {input.sorted_pixel_tsv} \
            --output {output.png_highlight} \
            --color_table $tmp_cmap \
            --plot_um_per_pixel {params.pixel_resolution} \
            --channel_list $channel_list  \
            --pcut 0.05  \
            --spcut 0.2

        """
# + ficture plot_pixel_full 
# \--input ../../output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/data/WT_37RF/analysis/nF20.d_12/nF20.d_12.decode.prj_12.r_4_5.pixel.sorted.tsv
# .gz
#  --color_table /vast/projects/CDFC/rlyu/Projects/Mouse_Myeloma_2025/output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/celltype.rgb.tsv 
#  --output ../../output/expr/proseg2.0.5/20251117_foci_manual_cutoff/ficture_out/predefine_model_m/data/WT_37RF/analysis/nF20.d_12/figure/nF20.d_12.decode.prj_12.r_4_5.pixel.png --plot_um_per_pixel 0.5 --full
