#!/bin/bash

#SBATCH -o /opt/mesh/eigg/sanket/cmdstan/Treg_dynamics/slurm_out/%j.%N.out
#SBATCH --error=/opt/mesh/eigg/sanket/cmdstan/Treg_dynamics/slurm_out/%j.%N.err_out
#SBATCH --get-user-env
#SBATCH -J sanket_Tregs
#SBATCH -D /opt/mesh/eigg/sanket/cmdstan/Treg_dynamics
#SBATCH --partition=general
#SBATCH --exclude=eigg,raasay,tiree,taransay
#SBATCH --ntasks=33
#SBATCG -N 2

while getopts m: flag
do
    case "${flag}" in
        m) modelname=${OPTARG};;
    esac
done


echo "stan_models/${modelname} sample num_samples=500 num_warmup=300 random seed=5689 id=1 data file=data/Treg_Newdata_S33.Rdump output file=save_csv/${modelname}_c.csv";
mpirun -vvvv stan_models/${modelname} sample num_samples=500 num_warmup=300 random seed=5689 id=1 data file=data/Treg_Newdata_S33.Rdump output file=save_csv/${modelname}_c1.csv &
mpirun -vvvv stan_models/${modelname} sample num_samples=500 num_warmup=300 random seed=5689 id=2 data file=data/Treg_Newdata_S33.Rdump output file=save_csv/${modelname}_c2.csv &
mpirun -vvvv stan_models/${modelname} sample num_samples=500 num_warmup=300 random seed=5689 id=3 data file=data/Treg_Newdata_S33.Rdump output file=save_csv/${modelname}_c3.csv 
