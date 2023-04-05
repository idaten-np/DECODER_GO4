

#if test -f "Go4AutoSave.root"; then
#	rm Go4AutoSave.root
#fi

file=khala36_run02_Eu152
#file=khala36_run01_Co60
for i in {0001..0100}
do
	sh ./condor_excutable.sh . go4analysis ${file}_${i}
done
