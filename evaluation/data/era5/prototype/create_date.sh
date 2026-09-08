#! /bin/sh -l

export scripts=/scratch4/NCEPDEV/land/data/evaluation/ush
export EXEC_DIR=../sorc

sdate=20240901
edate=20241231

rm -rf date_input.txt

while [ $sdate -le $edate ]; do
  
echo $sdate  >> date_input.txt

echo $sdate

sdate=`$scripts/finddate.sh $sdate d+1`

done
