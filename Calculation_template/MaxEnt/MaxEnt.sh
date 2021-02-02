#!/bin/bash
#

################################################################################
#                                USER SETTINGS                                 #
################################################################################
GENMAT=/home/petocchif/1_GW_EDMFT/nobackup/1_production/a_LiTi2O4               # Replace last field with project name

G_model=" " #  -M 0.5,1"
W_model=" " #  -m "$SRC"/Wmodel.dat"
U_model=" " #  -m ../../2_Wt2g/t2_new/Wt2g.dat_dos.dat"  #  -m "$SRC"/Wmodel.dat"




################################################################################
#                              DEFAULT ARGUMENTS                               #
################################################################################
err="0.008"                    # -e
mesh="6000"                    # -w
width="200"                    # -W
#FIELD=    MUST BE PROVIDED    # -f Available: "G" "W" "curlyU" "C" "M" "P"
#ORB=      MUST BE PROVIDED    # -o
JOR="1"                        # -j
SPIN="1"                       # -s Available: "1" "2"
#SOURCE=   MUST BE PROVIDED    # -i Available: "lat", "imp", "qmc_El"




################################################################################
#                             PROVIDED ARGUMENTS                               #
################################################################################
while getopts ":e:w:W:f:o:s:i:m:p:" o; do
   case ${o} in
      e)
         #echo ${OPTARG}
         err="${OPTARG}"
         ;;
      w)
         #echo ${OPTARG}
         mesh="${OPTARG}"
         ;;
      W)
         #echo ${OPTARG}
         width="${OPTARG}"
         ;;
      f)
         #echo ${OPTARG}
         FIELD="${OPTARG}"
         if [ "$FIELD"  == "G" ];      then RUNOPTIONS=" -s      "$err" -w "$mesh" -W "$width$G_model ; fi
         if [ "$FIELD"  == "W" ];      then RUNOPTIONS=" -S X -s "$err" -w "$mesh" -W "$width$W_model ; fi
         if [ "$FIELD"  == "curlyU" ]; then RUNOPTIONS=" -S X -s "$err" -w "$mesh" -W "$width$U_model ; fi
         if [ "$FIELD"  == "C" ];      then RUNOPTIONS=" -S X -s "$err" -w "$mesh" -W "$width$W_model ; fi
         if [ "$FIELD"  == "M" ];      then RUNOPTIONS=" -S X -s "$err" -w "$mesh" -W "$width$W_model ; fi
         if [ "$FIELD"  == "P" ];      then RUNOPTIONS=" -S X -s "$err" -w "$mesh" -W "$width$W_model ; fi
         if [[ "$FIELD"  = *PiJ* ]];   then RUNOPTIONS=" -S X -s "$err" -w "$mesh" -W "$width$U_model ; fi
         ;;
      o)
         #echo ${OPTARG}
         ORB="${OPTARG}"
         JOR="${OPTARG}"
         #if [ "$ORB"  != "t2g" ] && [ "$ORB"  != "eg" ]; then echo "Option Error - o" ; exit 1 ; fi
         ;;
      j)
         #echo ${OPTARG}
         JOR="${OPTARG}"
         ;;
      s)
         #echo ${OPTARG}
         SPIN="${OPTARG}"
         if [ "$SPIN"  != "1" ] && [ "$SPIN"  != "2" ]; then echo "Option Error - s" ; exit 1 ; fi
         ;;
      i)
         #echo ${OPTARG}
         SOURCE="${OPTARG}"
         #if [ "$SOURCE"  != "imp" ] && [ "$SOURCE"  != "lat" ]; then echo "Option Error - i" ; exit 1 ; fi
         ;;
      \? )
         echo "Invalid option: $OPTARG" 1>&2
         exit 1
         ;;
      : )
         echo "Invalid option: $OPTARG requires an argument" 1>&2
         exit 1
         ;;
    esac
done




################################################################################
#                               DATA EXTRACTION                                #
################################################################################
if [ "$FIELD"  == "G" ]; then
   #
   WORKDIR="MaxEnt_"${FIELD}${SOURCE}"_o"${ORB}"_s"${SPIN}
   mkdir $WORKDIR
   cd $WORKDIR
   JOBNAME=${FIELD}${ORB}"s"${SPIN}
   #
   DATA=${FIELD}${SOURCE}"_t_o"${ORB}"_s"${SPIN}".DAT"
   DATA_=$DATA
   cp ../${DATA} ./${DATA_}
   #
elif [ "$FIELD"  == "M" ]; then
   #
   WORKDIR="MaxEnt_"$FIELD$SOURCE
   mkdir $WORKDIR
   cd $WORKDIR
   JOBNAME=$FIELD$ORB$JOR
   #
   DATA=${FIELD}${SOURCE}"_w.DAT"
   DATA_=$DATA
   cp ../${DATA} ./${DATA_}
   #
else
   #
   WORKDIR="MaxEnt_"$FIELD$SOURCE"_o"$ORB$JOR
   mkdir $WORKDIR
   cd $WORKDIR
   JOBNAME=$FIELD$ORB$JOR
   #
   DATA=${FIELD}${SOURCE}"_w_("${ORB}","${JOR}")("${JOR}","${ORB}").DAT"
   DATA_=${FIELD}${SOURCE}"_w_"${ORB}${JOR}${JOR}${ORB}".DAT"
   cp ../${DATA} ./${DATA_}
   #
fi




################################################################################
#                                  PRINT INFOS                                 #
################################################################################
BIN=${GENMAT}/MaxEnt
echo
echo "Binary from: " ${BIN}
echo "Run options: "${RUNOPTIONS}
echo "File: "${DATA_}
echo "Jobname: "${JOBNAME}
echo




################################################################################
#                                 SCRIPT CREATION                              #
################################################################################
cat << EOF > submit_MaxEnt
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N   $JOBNAME
#$ -e   error.out
#$ -o   log.out
#$ -pe  smp 1
#$ -q   big.q

echo \$RUNOPTIONS
export PYTHONPATH=\${PYTHONPATH}:${BIN}/docopt/
export OMP_NUM_THREADS=1

mpiexec -np 1 python3.6  $BIN/bryan.py $RUNOPTIONS $DATA_ > job.out

EOF

qsub submit_MaxEnt
