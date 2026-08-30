#!/usr/bin/env bash
# Production template for g <= 2000 and all n >= 0.
# Run from the extracted high_genus_modular_checker directory.
set -euo pipefail

G=${G:-2000}
JOBS=${JOBS:-8}
TAIL_WORKERS=${TAIL_WORKERS:-32}
PRIMES=${PRIMES:-"5003 5009 5011 5021 5023 5039 5051 5059"}

mkdir -p tables logs

echo "[1/5] Building C++ programs"
g++ -O3 -std=c++17 -march=native \
  modular_low_truncation_table_large.cpp \
  -o modular_low_truncation_table_large
g++ -O3 -std=c++17 -march=native \
  check_modular_cover.cpp \
  -o check_modular_cover

echo "[2/5] Preparing W=14 table jobs"
: > jobs_G${G}_W14.txt
for p in ${PRIMES}; do
  s=$((p-1))
  table="tables/G${G}_W14_p${p}.bin"
  out="logs/G${G}_W14_p${p}.out"
  err="logs/G${G}_W14_p${p}.err"
  if [[ -s "${table}" ]] && grep -q '^PASS checksum=' "${out}" 2>/dev/null; then
    echo "reusing ${table}"
    continue
  fi
  cat >> jobs_G${G}_W14.txt <<EOF
./modular_low_truncation_table_large ${G} 14 ${s} ${p} ${table} >${out} 2>${err}
EOF
done

if [[ -s jobs_G${G}_W14.txt ]]; then
  if command -v parallel >/dev/null 2>&1; then
    echo "[3/5] Generating missing tables with GNU parallel, -j ${JOBS}"
    parallel --halt soon,fail=1 -j "${JOBS}" --joblog logs/table_jobs.tsv < jobs_G${G}_W14.txt
  else
    echo "GNU parallel is not installed; running the missing jobs sequentially"
    bash jobs_G${G}_W14.txt
  fi
else
  echo "[3/5] All requested tables already exist and have PASS logs"
fi

pass_count=0
prime_count=0
for p in ${PRIMES}; do
  prime_count=$((prime_count+1))
  if grep -q '^PASS checksum=' "logs/G${G}_W14_p${p}.out" 2>/dev/null; then
    pass_count=$((pass_count+1))
  fi
done
if [[ "${pass_count}" -ne "${prime_count}" ]]; then
  echo "ERROR: expected ${prime_count} PASS checksums, found ${pass_count}" >&2
  exit 1
fi

echo "[4/5] Finite ranges. Exit code 1 means unresolved pairs were written;"
echo "      this is expected while known zeros are still present."
set +e
./check_modular_cover "${G}" 10 0 49999 unresolved_a.tsv tables/G${G}_W14_p*.bin
status_a=$?
./check_modular_cover "${G}" 14 0 59999 unresolved_b.tsv tables/G${G}_W14_p*.bin
status_b=$?
set -e
if [[ ${status_a} -gt 1 || ${status_b} -gt 1 ]]; then
  echo "ERROR: a finite checker failed for a reason other than unresolved pairs" >&2
  exit 1
fi

echo "[5/5] Universal D1 tails"
python3 check_tail_prime_windows.py \
  --G "${G}" --Gamma 10 --tail-start 50000 --scan-end 2000000 \
  --infinite-start 2000000 --workers "${TAIL_WORKERS}" \
  | tee tail_a.log
python3 check_tail_prime_windows.py \
  --G "${G}" --Gamma 14 --tail-start 60000 --scan-end 4000000 \
  --infinite-start 4000000 --workers "${TAIL_WORKERS}" \
  | tee tail_b.log

cat <<EOF
DONE.
Now compare unresolved_a.tsv and unresolved_b.tsv with your exact known-zero lists.
If unexpected unresolved pairs remain, add more primes to PRIMES and rerun the
finite-check stage (existing table files can be reused).
EOF
