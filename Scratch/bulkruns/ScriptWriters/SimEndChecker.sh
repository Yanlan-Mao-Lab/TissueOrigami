n=180
i=140
while [ $i -lt $n ]
do
	printf -v CurrNo "%05d" $i
	echo $CurrNo
	tail -1 ../Run$CurrNo/Save_Frame	
	#sed -i 's/TimeStep(sec): 120/TimeStep(sec): 300/' ../ModelInputs/modelinput$CurrNo
	i=$(( $i + 1 ))
done
