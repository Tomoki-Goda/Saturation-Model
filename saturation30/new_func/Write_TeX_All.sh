

#alldir="ThetaOff"
#alldir="FirstRun"
alldir="NewRun"

gcc ./Utilities/Write_TeX.c -o Write_TeX -lm

for dir in  "BGK" #"GBS" "GBWIntegrated" "GBSPert" "GBSrfix" "GBSPertrfix" 
do

	./Write_TeX ./${alldir}/${dir}/Mass0.0-Qup10-*/result.txt \
		./${alldir}/${dir}/Mass0.0-Qup50-*/result.txt \
		./${alldir}/${dir}/Mass0.0-Qup100-*/result.txt \
		./${alldir}/${dir}/Mass0.0-Qup500-*/result.txt \
		./${alldir}/${dir}/massless.txt
		
	./Write_TeX ./${alldir}/${dir}/Mass0.0196-Qup10-*/result.txt \
		./${alldir}/${dir}/Mass0.0196-Qup50-*/result.txt \
		./${alldir}/${dir}/Mass0.0196-Qup100-*/result.txt \
		./${alldir}/${dir}/Mass0.0196-Qup500-*/result.txt \
		./${alldir}/${dir}/massive.txt
done

rm Write_TeX 
