#! /bin/bash

output_dir=output
bin_dir=../build/bin
data_dir=../data/meshes
camera_dir=../data/cameras

model_list=(\
 "killeroo_small_tri_clean"\
)

for model in ${model_list[@]}
do
  mkdir -p ${output_dir}/fig-perspective-distortion/${model}/low_fov
  ${bin_dir}/generate_perspective_distortion_figure \
    -i ${data_dir}/${model}_conf_simplified_with_uv.obj \
    -o ${output_dir}/fig-perspective-distortion/${model}/low_fov \
    --translation 0 --perspective_fov 25 --orthographic_fov 45

  mkdir -p ${output_dir}/fig-perspective-distortion/${model}/high_fov
  ${bin_dir}/generate_perspective_distortion_figure \
    -i ${data_dir}/${model}_conf_simplified_with_uv.obj \
    -o ${output_dir}/fig-perspective-distortion/${model}/high_fov \
    --translation -2.75 --perspective_fov 100 --orthographic_fov 20
done
