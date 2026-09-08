#! /bin/sh

export scripts=/scratch4/NCEPDEV/land/data/evaluation/ush
export EXEC_DIR=../sorc

sdate=20260101
edate=20260531

rm date_input.txt

while [ $sdate -le $edate ]; do
  
echo $sdate  >> date_input.txt

echo $sdate

sdate=`$scripts/finddate.sh $sdate d+1`

done
