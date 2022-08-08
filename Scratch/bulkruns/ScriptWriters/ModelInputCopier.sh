n=76
i=61000
#offset=16
#while [ $i -le $n ]
#do
	#printf -v CurrNo "%05d" $i
	#printf -v NewNo "%05d" $(( $i + $offset ))
	#cp  ../ModelInputs/modelinput$CurrNo ../ModelInputs/modelinput$NewNo
        #echo  ../ModelInputs/modelinput$CurrNo ../ModelInputs/modelinput$NewNo

	#i=$(( $i + 1 ))
#done

n=21
i=14
offset=1
while [ $i -le $n ]
do
	printf -v CurrNo "%05d" $i
	printf -v NewNo "%05d" $(( $i + $offset ))
	cp  ../ModelInputs/modelinput$CurrNo ../ModelInputs/modelinput$NewNo
        echo  ../ModelInputs/modelinput$CurrNo ../ModelInputs/modelinput$NewNo

	i=$(( $i + 1 ))
done

