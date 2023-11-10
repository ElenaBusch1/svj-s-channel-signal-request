cd 100xxx
for d in {100000..100047}; do
	echo "mv $d $[$d+415479]"
	mv $d $[$d+415479]
done
for d in {100048..100057}; do
	echo "mv $d $[$d+409914]"
	mv $d $[$d+409914]
done
cd ..

#for d in {515479..515526}; do
#	echo "mv $d $[$d-415479]"
#	mv $d $[$d-415479]
#done
#for d in {509962..509971}; do
#	echo "mv $d $[$d-409914]"
#	mv $d $[$d-409914]
#done

