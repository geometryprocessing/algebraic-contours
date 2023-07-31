#! /bin/bash

output_dir=output
bin_dir=../build/bin
data_dir=../data/meshes
camera_dir=../data/cameras

num_tests=26
mkdir -p ${output_dir}/fig-timing

model_list=($(ls ${data_dir}))

${bin_dir}/add_header_to_csv ${output_dir}/fig-timing

for file in ${model_list[@]}
do
  if [[ "$file" = *".obj" ]]; then
    model=${file%.*}
    ${bin_dir}/generate_timing_data \
      -i ${data_dir}/${model}.obj \
      -o ${output_dir}/fig-timing \
      --num_tests ${num_tests}
  fi
done