# Code for ShapeMotifEM benchmarking

For each TF analyzed we ran the following basic example code on Expanse
to compare AUPRs arrived at by motif models trained by ShapeME and ShapeMotifEM
using identical inputs.

```bash
singularity exec -B ${data_dir}:${data_dir} \
    /home/jwschroeder/appliances/ShapeMotifEM/shapemotifem.sif \
    Rscript /src/shapemotifem_analysis_code/learn_motifs.R \
    -f {train_seqs} \
    -d ${data_dir}

singularity exec -B ${data_dir}:${data_dir} \
    /home/jwschroeder/appliances/ShapeMotifEM/shapemotifem.sif \
    Rscript /src/shapemotifem_analysis_code/standardize_shapes.R \
    --fa_file ${test_seqs} \
    --data_direc ${data_dir} \
    --score_file ${test_scores} \
    --motifs_file ${shapemotifem_results}

singularity exec -B ${data_dir}:${data_dir} \
    /home/jwschroeder/appliances/ShapeMotifEM/shapemotifem.sif \
    Rscript /src/shapemotifem_analysis_code/get_distances.R \
    --data_direc ${data_dir} \
    --prefix "test_${dataset_prefix}" \
    --std_data ${standardize_shapes_result}
```
