#Plot piecharts from LEA results as maps
#from the results in 'Structure' folder of SambaR results
#Tomas Fer, 2025-2026
#tomas.fer@natur.cuni.cz

#REQUIREMENTS
#GhostScript (gs)
#R + packages: maps, plotrix, RColorBrewer
#pdfbox (automatically downloaded from Apache)

#takes output of LEAqmatrix for K=2 to k (defined below)
#adds localities column to the LEA output (requires localitiesCoor.txt with 3 columns: locID (the same as in qmatrix), x, y)
#uses plotLEA_maps.R to plot the piecharts on the map (works well with a single data point per locality)
#combines all PDFs to a single file using 'pdfbox'

#define max K
k=${1}
for i in $(seq 2 $k); do
	echo "Plotting K${i}"
	#join LEAqmatrix output file with coordinate file (both have sample name as a first column)
	join <(sort localitiesCoor.txt) <(sort Structureplot.LEAqmatrix_K${i}.txt) | tr ' ' '\t' > K${i}.txt
	#create a header
	echo -ne "name\tx\ty\t" > header
	for j in $(seq 1 $i); do echo -ne "c${j}\t" >> header; done
	echo >> header
	cat header K${i}.txt > tmp && mv tmp K${i}.txt
	rm header
	#run map plotting script
	./plotLEA_maps.R ${i}
	#crop resulting PDF (keep margin 10 point)
	gs -o null -sDEVICE=bbox K${i}.pdf 2>out > gs.log #reports crop box coordinates
	echo >> gs.log
	#take the four numbers and add/subtract a value ('add')
	add=10
	crop1=$(grep HiRes out | cut -d' ' -f2) #left margin
	crop1x=$(echo "$crop1 - $add" | bc)
	crop2=$(grep HiRes out | cut -d' ' -f3) #bottom margin
	crop2x=$(echo "$crop2 - $add" | bc)
	crop3=$(grep HiRes out | cut -d' ' -f4) #right margin
	crop3x=$(echo "$crop3 + $add" | bc)
	crop4=$(grep HiRes out | cut -d' ' -f5) #upper margin
	crop4x=$(echo "$crop4 + $add" | bc)
	crop=$(echo $crop1x $crop2x $crop3x $crop4x)
	gs -o K${i}crop.pdf -sDEVICE=pdfwrite -dAutoRotatePages=/None -dUseCropBox=true -c "[/CropBox [$crop] /PAGES pdfmark" -f K${i}.pdf >> gs.log
	echo >> gs.log
	mv K${i}crop.pdf K${i}.pdf
	rm out
done
#merge all PDFs into a single file
echo "Merging PDFs"
#Check the newest PDFbox version and silently download it
pdfboxver=$(wget -q -O- https://downloads.apache.org/pdfbox/ | grep "2\.0\." | cut -d'"' -f6 | sed 's/.$//')
#download the newest version and rename
wget -q https://downloads.apache.org/pdfbox/${pdfboxver}/pdfbox-app-${pdfboxver}.jar
mv pdfbox-app-${pdfboxver}.jar pdfbox.jar
#merge K*.pdf to a single PDF
java -jar pdfbox.jar PDFMerger K*.pdf ${2}_K2-${k}.pdf

