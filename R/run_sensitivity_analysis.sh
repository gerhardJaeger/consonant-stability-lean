for i in $(seq 1 4); do
  for j in $(seq 1 12); do
    echo "Running i=$i j=$j"
    Rscript OU_sensitivity_analysis.r "$i" "$j"
  done
done