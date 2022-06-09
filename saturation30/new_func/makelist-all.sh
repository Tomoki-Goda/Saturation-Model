

#alldir="ThetaOff"
#alldir="FirstRun"
alldir="NewRun"

gcc makelist.c -o makelist -lm

for dir in "GBS" "GBWIntegrated" "GBSPert" "GBSrfix" "GBSPertrfix" 
do

	./makelist ./${alldir}/${dir}/Mass0.0-Qup10-*/result.txt \
		./${alldir}/${dir}/Mass0.0-Qup50-*/result.txt \
		./${alldir}/${dir}/Mass0.0-Qup100-*/result.txt \
		./${alldir}/${dir}/Mass0.0-Qup500-*/result.txt \
		./${alldir}/${dir}/massless.txt
		
	./makelist ./${alldir}/${dir}/Mass0.0196-Qup10-*/result.txt \
		./${alldir}/${dir}/Mass0.0196-Qup50-*/result.txt \
		./${alldir}/${dir}/Mass0.0196-Qup100-*/result.txt \
		./${alldir}/${dir}/Mass0.0196-Qup500-*/result.txt \
		./${alldir}/${dir}/massive.txt
done
