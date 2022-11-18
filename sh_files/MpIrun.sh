#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/slurm_out/%j.%N.out
#SBATCH --error=/opt/mesh/eigg/sanket/slurm_out/%j.%N.err_out
#SBATCH --get-user-env
#SBATCH -J test_sanket
#SBATCH -D /opt/mesh/eigg/sanket/cmdstan/Treg_dynamics
#SBATCH --nodelist=tiree,raasay
#SBATCH --ntasks=33

while getopts m: flag
do
    case "${flag}" in
        m) modelname=${OPTARG};;
    esac
done


echo "stan_models/MAP_${modelname}_naiTreg sample num_samples=500 num_warmup=300 data file=data/Treg_data.Rdump output file=save_csv/${modelname}_s33.csv";
mpirun -vvvv stan_models/MAP_${modelname}_naiTreg sample num_samples=500 num_warmup=300 data file=data/Treg_data.Rdump output file=save_csv/${modelname}_s33.csv
