#!/bin/bash
#
QUEUES="new.q"
#
################################################################################
#                                USER SETTINGS                                 #
################################################################################
G_model=" " #  -M 0.5,1"
W_model=" " #  -m "$SRC"/Wmodel.dat"
U_model=" " #  -m ../../2_Wt2g/t2_new/Wt2g.dat_dos.dat"  #  -m "$SRC"/Wmodel.dat"




################################################################################
#                              DEFAULT ARGUMENTS                               #
################################################################################
err="0.008"                    # -e
mesh="6000"                    # -w
width="200"                    # -W
#FIELD=    MUST BE PROVIDED    # -f Available: "G" "W" "curlyU" "C" "M" "P" "S"(only for SOURCE=qmc_El)
#ORB=      MUST BE PROVIDED    # -o
JOR="1"                        # -j
SPIN="1"                       # -s Available: "1" "2"
#SOURCE=   MUST BE PROVIDED    # -i Available: "lat", "imp", "qmc_El"
BINPATH=/home/petocchif/1_GW_EDMFT/GenericMaterial/Calculation_template




################################################################################
#                             PROVIDED ARGUMENTS                               #
################################################################################
while getopts ":e:w:W:f:o:j:s:i:p:" o; do
   case ${o} in
      e) #ERROR
         #echo ${OPTARG}
         err="${OPTARG}"
         ;;
      w) #NUMBER OF REAL FREQUENCY POINTS
         #echo ${OPTARG}
         mesh="${OPTARG}"
         ;;
      W) #WIDTH OF THE REAL FREQUENCY MESH
         #echo ${OPTARG}
         width="${OPTARG}"
         ;;
      f) #INPUT FIELD
         #echo ${OPTARG}
         FIELD="${OPTARG}"
         if [ "$FIELD"  == "G" ];      then RUNOPTIONS=" -S F -s "$err" -w "$mesh" -W "$width$G_model ; fi
         if [ "$FIELD"  == "S" ];      then RUNOPTIONS=" -S F -s "$err" -w "$mesh" -W "$width$G_model ; fi
         if [ "$FIELD"  == "W" ];      then RUNOPTIONS=" -S X -s "$err" -w "$mesh" -W "$width$W_model ; fi
         if [ "$FIELD"  == "curlyU" ]; then RUNOPTIONS=" -S X -s "$err" -w "$mesh" -W "$width$U_model ; fi
         if [ "$FIELD"  == "C" ];      then RUNOPTIONS=" -S X -s "$err" -w "$mesh" -W "$width$W_model ; fi
         if [ "$FIELD"  == "M" ];      then RUNOPTIONS=" -S X -s "$err" -w "$mesh" -W "$width$W_model ; fi
         if [ "$FIELD"  == "P" ];      then RUNOPTIONS=" -S X -s "$err" -w "$mesh" -W "$width$W_model ; fi
         if [[ "$FIELD"  = *PiJ* ]];   then RUNOPTIONS=" -S X -s "$err" -w "$mesh" -W "$width$U_model ; fi
         ;;
      o) #DIAGONAL ORBITAL INDEX
         #echo ${OPTARG}
         ORB="${OPTARG}"
         JOR="${OPTARG}"
         #if [ "$ORB"  != "t2g" ] && [ "$ORB"  != "eg" ]; then echo "Option Error - o" ; exit 1 ; fi
         ;;
      j) #OFF-DIAGONAL ORBITAL INDEX
         #echo ${OPTARG}
         JOR="${OPTARG}"
         ;;
      s) #SPIN INDEX
         #echo ${OPTARG}
         SPIN="${OPTARG}"
         if [ "$SPIN"  != "1" ] && [ "$SPIN"  != "2" ]; then echo "Option Error - s" ; exit 1 ; fi
         ;;
      i) #INPUT FIELD SUFFIX
         #echo ${OPTARG}
         SOURCE="${OPTARG}"
         ;;
      p) #PATH TO MAXENT FOLDER
         #echo ${OPTARG}
         BINPATH="${OPTARG}"
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
if [ "$FIELD"  == "G" ] || [ "$FIELD"  == "S" ]; then
   #
   WORKDIR="MaxEnt_"${FIELD}${SOURCE}"_o"${ORB}"_s"${SPIN}
   mkdir $WORKDIR
   cd $WORKDIR
   JOBNAME=${FIELD}${SOURCE}"_"${ORB}"s"${SPIN}
   #
   DATA=${FIELD}${SOURCE}"_t_o"${ORB}"_s"${SPIN}".DAT"
   DATA_=$DATA
   cp ../${FIELD}${SOURCE}/${DATA} ./${DATA_}
   #
elif [ "$FIELD"  == "M" ]; then
   #
   WORKDIR="MaxEnt_"${FIELD}${SOURCE}
   mkdir $WORKDIR
   cd $WORKDIR
   JOBNAME=${FIELD}${SOURCE}"_"${ORB}${JOR}
   #
   DATA=${FIELD}${SOURCE}"_w.DAT"
   DATA_=$DATA
   cp ../${FIELD}${SOURCE}/${DATA} ./${DATA_}
   #
else
   #
   WORKDIR="MaxEnt_"${FIELD}${SOURCE}"_o"${ORB}${JOR}
   mkdir $WORKDIR
   cd $WORKDIR
   JOBNAME=${FIELD}${SOURCE}"_"${ORB}${JOR}
   #
   DATA=${FIELD}${SOURCE}"_w_("${ORB}","${JOR}")("${JOR}","${ORB}").DAT"
   DATA_=${FIELD}${SOURCE}"_w_"${ORB}${JOR}${JOR}${ORB}".DAT"
   cp ../${FIELD}${SOURCE}/${DATA} ./${DATA_}
   #
fi




################################################################################
#                                  PRINT INFOS                                 #
################################################################################
BIN=${BINPATH}/MaxEnt
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
#SBATCH --job-name=${JOBNAME}
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --account=unifr
#SBATCH --partition=new

echo \$RUNOPTIONS
export PYTHONPATH=\${PYTHONPATH}:${BIN}/docopt/
export OMP_NUM_THREADS=1

srun python3.6  ${BIN}/bryan.py ${RUNOPTIONS} ${DATA_} > job.out 2> err.out

EOF

sbatch submit_MaxEnt
