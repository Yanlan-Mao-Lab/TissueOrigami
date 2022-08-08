n=1
i=162

#109 -> 124
offset=1
while [ $i -le $n ]
do
	printf -v CurrNo "%05d" $i
	printf -v NewNo "%05d" $(( $i + $offset ))
        echo  ../Run$CurrNo/ ../Run$NewNo/

	cp  ../Run$CurrNo/Save* ../Run$NewNo/
        cp  ../Run$CurrNo/Out ../Run$NewNo/

	i=$(( $i + 1 ))
done


