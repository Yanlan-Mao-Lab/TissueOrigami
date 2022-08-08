RunNo=7000

while [ $RunNo -le 7004 ]
do
	printf -v CurrNo "%05d" $RunNo	
	#qsub -cwd -l h_rt=5:00:00 -l memory=16G -l tmpfs=6G -pe smp 4 -S /bin/bash -N TF$CurrNo  -wd /home/ucgamto/Scratch/Projects/TissueFolding/bulkruns/SubmissionFolder/ /home/ucgamto/Scratch/Projects/TissueFolding/bulkruns/Run$CurrNo/Job$CurrNo	
	qsub -cwd -l h_rt=12:00:00 -l memory=16G -l tmpfs=6G -pe smp 4 -S /bin/bash -N TF$CurrNo  -wd /home/ucbpnkh/Scratch/bulkruns/SubmissionFolder/ /home/ucbpnkh/Scratch/bulkruns/Run$CurrNo/Job$CurrNo 
RunNo=$(( $RunNo + 1 ))
done
