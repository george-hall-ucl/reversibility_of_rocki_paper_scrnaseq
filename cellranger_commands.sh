# (C) University College London 2020-2022. This software is licenced under the
# terms of the GNU GENERAL PUBLIC LICENSE Version 3. See COPYING.txt for the
# licence details.

# We used cellranger v5.0.1
# "SEQ_ID" and "RUN_ID" are used as placeholders for the sake of neatness.

# effect-of-rocki-three-samples.csv:
#
# Lane,Sample,Index
# *,1CTRL6D,SI-TT-A1
# *,1ROCKi6D,SI-TT-B1
# *,1CTRL12D,SI-TT-C1
# *,1ROCKi12D,SI-TT-D1
# *,2CTRL6D,SI-TT-E1
# *,2ROCKi6D,SI-TT-F1
# *,2CTRL12D,SI-TT-G1
# *,2ROCKi12D,SI-TT-H1
# *,3CTRL6D,SI-TT-A2
# *,3ROCKi6D,SI-TT-B2
# *,3CTRL12D,SI-TT-C2
# *,3ROCKi12D,SI-TT-D2

cellranger mkfastq --run SEQ_ID --csv=effect-of-rocki-three-samples.csv
# produces the directory "RUN_ID"

# download reference if necessary:
# wget https://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
# tar -xzvf refdata-cellranger-GRCh38-3.0.0.tar.gz

for SAMPLE in "1CTRL6D" "1ROCKi6D" "1CTRL12D" "1ROCKi12D" "2CTRL6D" "2ROCKi6D" "2CTRL12D" "2ROCKi12D" "3CTRL6D" "3ROCKi6D" "3CTRL12D" "3ROCKi12D"; do
    cellranger count --id=$SAMPLE --fastqs=RUN_ID/outs/fastq_path --sample=$SAMPLE --transcriptome=refdata-cellranger-GRCh38-3.0.0
done

# aggregation_info.csv:
#
# library_id,molecule_h5,control_or_treated,timepoint,treatment_and_timepoint
# 1CTRL6D,1CTRL6D/outs/molecule_info.h5,control,6D,CTRL6D
# 1ROCKi6D,1ROCKi6D/outs/molecule_info.h5,treated,6D,ROCKi6D
# 1CTRL12D,1CTRL12D/outs/molecule_info.h5,control,12D,CTRL12D
# 1ROCKi12D,1ROCKi12D/outs/molecule_info.h5,treated,12D,ROCKi12D
# 2CTRL6D,2CTRL6D/outs/molecule_info.h5,control,6D,CTRL6D
# 2ROCKi6D,2ROCKi6D/outs/molecule_info.h5,treated,6D,ROCKi6D
# 2CTRL12D,2CTRL12D/outs/molecule_info.h5,control,12D,CTRL12D
# 2ROCKi12D,2ROCKi12D/outs/molecule_info.h5,treated,12D,ROCKi12D
# 3CTRL6D,3CTRL6D/outs/molecule_info.h5,control,6D,CTRL6D
# 3ROCKi6D,3ROCKi6D/outs/molecule_info.h5,treated,6D,ROCKi6D
# 3CTRL12D,3CTRL12D/outs/molecule_info.h5,control,12D,CTRL12D
# 3ROCKi12D,3ROCKi12D/outs/molecule_info.h5,treated,12D,ROCKi12D

cellranger aggr --id=effect_of_rocki_3_donors --csv=aggregation_info.csv --normalize=none
