for f in 1000*; do
	echo "cp condor/$f/log.generate $f/log.generate"
	cp condor/$f/log.generate $f/log.generate
done
