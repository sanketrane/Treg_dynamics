
while getopts m: flag
do
    case "${flag}" in
        m) modelname=${OPTARG};;
    esac
done

for i in 1 2 3 4
    do
      ./models/${modelname} sample num_warmup=500 num_samples=500 data file=datafiles/MZ_data.Rdump \
      output file=save_csv/${modelname}_${i}.csv &
    done
