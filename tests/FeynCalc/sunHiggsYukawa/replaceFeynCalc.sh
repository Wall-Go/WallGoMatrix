#!/bin/bash

# Check if a file was provided
if [ $# -ne 1 ]; then
  echo "Usage: $0 <input_file>"
  exit 1
fi

input_file="$1"

# Extract filename parts
filename=$(basename -- "$input_file")
dirname=$(dirname -- "$input_file")
extension="${filename##*.}"
basename="${filename%.*}"

# Create new output filename
newExtension="_replaced"
# new_extension = ""
output_file="${dirname}/${basename}${newExtension}.${extension}"

# Perform replacements using Perl
perl -pe '
  s/IndexDelta\[\s*([^,\]]+)\s*,\s*([^\]]+)\s*\]/SDF[$1, $2]/g;
  s/SUNF\[\s*([^,\]]+)\s*,\s*([^,\]]+)\s*,\s*([^,\]]+)\s*,\s*([^\]]+)\s*\]/SUNF[$1, $3, Adjoint\$1] SUNF[Adjoint\$1, $4, $2]/g;
  s/SUNT\[\s*([^,\]]+)\s*,\s*([^,\]]+)\s*,\s*([^\]]+)\s*\]/SUNTF[{$1}, $2, $3]/g;
  s/IndexSum\[\s*(.*?),\s*\{Fundamental\$1,\s*1,\s*3\}\s*\]/$1/g;
  s/ScalarProduct/FAScalarProduct/g;
  s/MetricTensor/FAMetricTensor/g;
  s/FourVector/FAFourVector/g;
  s/NonCommutative/FANonCommutative/g;
  s/DiracSpinor/FADiracSpinor/g;
  s/PropagatorDenominator/FAPropagatorDenominator/g;
  s/PolarizationVector/FAPolarizationVector/g;
  s/ChiralityProjector/FAChiralityProjector/g;
  s/DiracSlash/FADiracSlash/g;
  s/DiracMatrix/FADiracMatrix/g;
  s/GaugeXi/FAGaugeXi/g;
' "$input_file" > "$output_file"

echo "Replacements complete. Output written to: $output_file"

# Copy associated .gen and .pars files
for ext in pars; do
  sidefile="${dirname}/${basename}.${ext}"
  if [ -f "$sidefile" ]; then
    cp "$sidefile" "${dirname}/${basename}${newExtension}.${ext}"
    echo "Copied $sidefile to ${basename}${newExtension}.${ext}"
  fi
done
