#!/bin/bash
set -euo pipefail
export LC_ALL=C

INPUT1="data/HGMTPointVac"
INPUT2="simulation_materials/kapton_effs.csv"
INPUT3="data/"
taus=(127 102 76 51)
alphas=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)

outdir="data"
mkdir -p "$outdir"
sanitize(){ echo "$1" | sed 's/\./p/g'; }

# preflight checks
[[ -f "$INPUT1" ]] || { echo "ERROR: INPUT1 not found: $INPUT1"; exit 1; }
[[ -f "$INPUT2" ]] || { echo "ERROR: INPUT2 not found: $INPUT2"; exit 1; }
[[ -d "$INPUT3" ]] || { echo "ERROR: OUTPUT DIR not found: $INPUT3"; exit 1; }
command -v make >/dev/null || { echo "ERROR: 'make' not found"; exit 1; }

for tau in "${taus[@]}"; do
  tau_tag=$(sanitize "$tau")
  outfile="${outdir}/tau_${tau_tag}.csv"
  echo "alpha,total_hits,total_scatters,metric" > "$outfile"

  for alpha in "${alphas[@]}"; do
    alpha_tag=$(sanitize "$alpha")
    echo "==> τ=${tau}, α=${alpha}"

    # Force rebuild with new macros (show output for debugging)
    make -B CFLAGS+=" -DTAU=${tau} -DALPHA=${alpha}" \
         CXXFLAGS+=" -DTAU=${tau} -DALPHA=${alpha}"

    # Confirm the executable exists
    [[ -x ./hgmt_lor_creator ]] || { echo "ERROR: build produced no ./hgmt_lor_creator"; exit 1; }

    log="log_tau${tau_tag}_alpha${alpha_tag}.txt"
    echo "[run] ./hgmt_lor_creator $INPUT1 $INPUT2 $INPUT3"
    # See output live AND save to log
    set +e
    ./hgmt_lor_creator "$INPUT1" "$INPUT2" "$INPUT3" 2>&1 | tee "$log"
    status=${PIPESTATUS[0]}
    set -e
    if [[ $status -ne 0 ]]; then
      echo "ERROR: hgmt_lor_creator exited with code $status (see $log)"
      exit $status
    fi

    # Parse totals (lines begin with 'total scatters:' / 'total hits:')
    scat=$(awk '/^total[[:space:]]+scatters:/ {print $3}' "$log" | tail -n1)
    hits=$(awk '/^total[[:space:]]+hits:/     {print $3}' "$log" | tail -n1)

    # Basic sanity
    if [[ -z "${scat}" || -z "${hits}" ]]; then
      echo "ERROR: Could not parse totals from $log"
      exit 1
    fi

    metric=$(awk -v h="$hits" -v s="$scat" 'BEGIN{ if(s>0) printf("%.6f", h/s); else print "NaN" }')
    echo "${alpha},${hits},${scat},${metric}" >> "$outfile"
  done
done

echo "Done. CSVs per τ in: $outdir/"
