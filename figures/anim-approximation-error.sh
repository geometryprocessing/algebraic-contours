#! /bin/bash

output_dir=output
bin_dir=../build/bin
data_dir=../data/meshes
camera_dir=../data/cameras

model_list=(\
  "fertility_tri_clean"\
)

for model in ${model_list[@]}
do
  mkdir -p ${output_dir}/anim-approximation-error/${model}
  ${bin_dir}/animate_approximation_error \
    -i ${data_dir}/${model}_conf_simplified_with_uv.obj \
    -o ${output_dir}/anim-approximation-error/${model}
done

for i in {0..360}
do
  compare ${output_dir}/projected_planar_frame_${i}.png \
    ${output_dir}/unprojected_planar_frame_${i}.png \
    ${output_dir}/difference_frame_${i}.png
done