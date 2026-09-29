#! /bin/sh -l

export scripts=/scratch2/NCEPDEV/land/data/evaluation/ush
export EXEC_DIR=../sorc

sdate=20150101
edate=20221231

#sdate=20150101
#edate=20150101

while [ $sdate -le $edate ]; do
   year=`echo $sdate |cut -c1-4`
   mon=`echo $sdate |cut -c5-6`
   day=`echo $sdate |cut -c7-8`
   julian=`$scripts/julian.sh $sdate`

   JDA=$julian

   if [ $julian -lt 10 ]
   then
      JDA=00$julian
   fi

   if [ $julian -ge 10 ] && [ $julian -lt 100 ]
   then
      JDA=0$julian
   fi

datestring=$year$JDA

\rm date_input.txt
  
echo $datestring > date_input.txt
echo $sdate >> date_input.txt
echo $year >> date_input.txt
echo $mon >> date_input.txt
echo $day >> date_input.txt

${EXEC_DIR}/regrid_ims_cgrid.exe
echo $sdate

sdate=`$scripts/finddate.sh $sdate d+1`

done
