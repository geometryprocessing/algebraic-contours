#! /bin/bash

output_dir=output
bin_dir=../build/bin
data_dir=../data/meshes
camera_dir=../data/cameras

model_list=(\
 "ogre_controlmesh-cleaned-3"\
)

for model in ${model_list[@]}
do
  for weight in 0.1 1 10 100
  do
    mkdir -p ${output_dir}/fig-smoothness/${model}/weight_${weight}
    ${bin_dir}/generate_smoothness_figure \
      -i ${data_dir}/${model}_conf_simplified_with_uv.obj \
      -o ${output_dir}/fig-smoothness/${model}/weight_${weight} \
      --weight ${weight}
  done
done
