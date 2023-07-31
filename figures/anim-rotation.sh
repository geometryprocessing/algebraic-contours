#! /bin/bash

output_dir=output
bin_dir=../build/bin
data_dir=../data/meshes
camera_dir=../data/cameras

model_list=(\
 "bigguy-cleaned-2"\
 "bob_quad_tri_clean"\
 "pipes_tri_clean"\
 "ribbon_closed_tri_clean"\
 "spot_quadrangulated_tri_clean"\
)

for model in ${model_list[@]}
do
  mkdir -p ${output_dir}/anim-rotation/${model}
  ${bin_dir}/animate_rotation \
    -i ${data_dir}/${model}_conf_simplified_with_uv.obj \
    -o ${output_dir}/anim-rotation/${model}

  file_list=($(ls ${output_dir}/anim-rotation/${model}))
  for file in ${file_list[@]}
  do
    if [[ "$file" = *".svg" ]]; then
      frame=${file%.*}
      convert -size 1600x1600 ${output_dir}/anim-rotation/${model}/${frame}.svg \
        ${output_dir}/anim-rotation/${model}/${frame}.png
    fi
  done
done
