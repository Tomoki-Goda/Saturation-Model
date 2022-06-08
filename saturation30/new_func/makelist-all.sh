

#alldir="ThetaOff"
alldir="FirstRun"

for dir in "GBS" "GBWIntegrated" "GBSPert" "GBSrfix" "GBSPertrfix" 
do

	./makelist ./${alldir}/${dir}/10masslessLCB/result.txt \
		./${alldir}/${dir}/50masslessLCB/result.txt \
		./${alldir}/${dir}/100masslessLCB/result.txt \
		./${alldir}/${dir}/500masslessLCB/result.txt \
		./${alldir}/${dir}/massless.txt
		
	./makelist ./${alldir}/${dir}/10massiveLCB/result.txt \
		./${alldir}/${dir}/50massiveLCB/result.txt \
		./${alldir}/${dir}/100massiveLCB/result.txt \
		./${alldir}/${dir}/500massiveLCB/result.txt \
		./${alldir}/${dir}/massive.txt
done
