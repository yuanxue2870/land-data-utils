#!/bin/bash

# A script to download all files from the NOAA G02158 unmasked directory.
# This script uses wget's recursive, no-parent, and reject features.

export scripts=/scratch4/NCEPDEV/land/data/evaluation/ush

folder=09_Sep
sdate=20250901
edate=20250901

while [ $sdate -le $edate ]; do

echo $sdate

wget --recursive -nH --cut-dirs=5 --no-parent --accept=.tar https://noaadata.apps.nsidc.org/NOAA/G02158/unmasked/2025/$folder/SNODAS_unmasked_$sdate.tar

tar -xvf SNODAS_unmasked_$sdate.tar
gunzip *.gz

for filename in zz_ssmv*$sdate*.dat; do
    # Skip if no match
    [[ -e "$filename" ]] || continue
    file="${filename%.*}"  # Remove extension
    echo "$file"
cat <<EOF > $file.hdr
ENVI
samples = 8192
lines = 4096
bands = 1
header offset = 0
file type = ENVI Standard
data type = 2
interleave = bsq
byte order = 1
EOF

done

sdate=`$scripts/finddate.sh $sdate d+1`

done

echo "Download process complete."
