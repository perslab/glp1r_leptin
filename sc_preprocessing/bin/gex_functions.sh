run_bcl2fastq() {
    local flowcell_dir=$1
    local output_dir=$2
    local sample_sheet=$3
    local tenx_type=$4
    local interop_dir=./stats

    local bcl_cmd="bcl2fastq --create-fastq-for-index-reads \
                              --minimum-trimmed-read-length=10 \
                              --mask-short-adapter-reads=10 \
                              --ignore-missing-positions \
                              --ignore-missing-controls \
                              --ignore-missing-filter \
                              --ignore-missing-bcls \
                              -r 10 -p 40 -w 10 \
                              -R ${flowcell_dir} \
                              --output-dir=${output_dir} \
                              --interop-dir=${interop_dir} \
                              --sample-sheet=${sample_sheet}"

    if [ "${tenx_type}" == "3.0" ]; then
        bcl_cmd+=" --barcode-mismatches 0"
    else
        bcl_cmd+=" --use-bases-mask=Y28,I10,I10,Y90"
    fi

    echo "Running bcl2fastq with command: $bcl_cmd"
    eval "$bcl_cmd"
}

run_cellranger() {
    local sample_id=$1
    local fastq_dir=$2
    local cores=$3
    local seq_type=$4
    local outpath=$5
    local transcriptome

    if [ "${seq_type}" == "tapseq" ]; then
        transcriptome="/projects/bernd/SCOP_2022_0188/mouse_optimized_v147"
    elif [ "${seq_type}" == "mouse" ]; then
        transcriptome="${outpath}/transcriptomes/mouse_mm10_optimized_v1/"
    else 
        transcriptome="${outpath}/transcriptomes/human_grch38_optimized_v1/"
    fi
    
    echo "Checking transcriptome path: ${transcriptome}"
    if [ $? -ne 0 ]; then
       echo "Failed to access transcriptome directory at ${transcriptome}."
        exit 1
    fi


    local cellranger_cmd="cellranger count --id=${sample_id} \
                        --transcriptome=${transcriptome} \
                        --fastqs=${fastq_dir} \
                        --sample=${sample_id}E \
                        --localcores=${cores}"

    local gex_counts_dir="${outpath}/gex_counts"

    # Check if the directory exists. If not, create it.
    if [ ! -d "${gex_counts_dir}" ]; then
        echo "Directory ${gex_counts_dir} does not exist. Creating now."
        mkdir -p "${gex_counts_dir}"
        if [ $? -ne 0 ]; then
            echo "Failed to create directory ${gex_counts_dir}."
            exit 1
        fi
    fi

    echo "Running Cell Ranger with command: ${cellranger_cmd}"
    cd "${gex_counts_dir}" || exit
    eval "${cellranger_cmd}"
    
    
    # New part: Move the raw_feature_bc.h5 file to the raw_h5_files directory
    local h5_output_dir="${outpath}/raw_h5_files"
    mkdir -p "${h5_output_dir}"
    local h5_file="${outpath}/gex_counts/${sample_id}/outs/raw_feature_bc_matrix.h5"
    local new_h5_file_name="${h5_output_dir}/${sample_id}_raw_feature_bc_matrix.h5"
    
    if [ -f "${h5_file}" ]; then
        cp "${h5_file}" "${new_h5_file_name}"
        echo "Exported h5 file to ${new_h5_file_name}"
    else
        echo "Error: h5 file does not exist at ${h5_file}"
    fi
}


run_kite() {
    local fastq_dir=$1
    local hto_dir=$2
    local feature_file=$3
    local cores=$4

    mkdir -p "${hto_dir}"
    echo "HTO MAPPING FOUND IN ${feature_file}"

    local fastqs=($(find "${fastq_dir}" -type f -name "*H*R*gz" | sort))
    if (( ${#fastqs[@]} > 0 )); then
        echo "HTO FASTQ FILES: ${fastqs[*]}"
        local kb_ref_cmd="kb ref -i ${hto_dir}/mismatch.idx -f1 ${hto_dir}/mismatch.fa -g ${hto_dir}/t2g.txt --workflow kite ${feature_file} --overwrite"
        local kb_count_cmd="kb count -i ${hto_dir}/mismatch.idx -g ${hto_dir}/t2g.txt -x 10xv3 --workflow kite -t ${cores} -o ${hto_dir} --overwrite --h5ad"

        echo "Running Kite with command: ${kb_ref_cmd}"
        eval "${kb_ref_cmd}"

        echo "Running Kite with command: ${kb_count_cmd}"
        eval "${kb_count_cmd}" "${fastqs[@]}"
    else
        echo "No HTO FASTQ files found, skipping."
    fi
}
