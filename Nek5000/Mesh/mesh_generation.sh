#!/bin/bash

outdir="IR21"
rm -rf ${outdir}
mkdir ${outdir}

############
# GEOMETRY #
############

# INPUTS
geo_json="./obstacle_geo/sample_io/geom_C4.json"
geo_case="bl" # or channel

# OUTPUTS
geo_out="${outdir}/IR21.geo"
mesh_info_json="${outdir}/mesh_info.json"
bc_json="${outdir}/bc.json"

# Create obstacle geometry
python3 obstacle_geo/obstacle_tools.py \
    ${geo_json} \
    ${geo_case} \
    ${geo_out} \
    ${mesh_info_json} \
    ${bc_json}

