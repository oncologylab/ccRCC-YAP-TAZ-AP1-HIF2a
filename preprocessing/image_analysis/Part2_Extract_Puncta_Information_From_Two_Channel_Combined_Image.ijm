// Puncta quantification using analyze particles optimized for SORA 60x Confocal images using Nikon Spinning Disk
run("Close All");
run("Set Measurements...", "area mean standard min center limit display redirect=None decimal=3");
Input = getDirectory("Choose a Directory "); //Select the Combined Images folder from the experiment of interest outputted in Part2
list = getFileList(Input);
Array.print(list);
MOD1 = 20;
MOD2 = 35;

tmp = Input + "Channel1_Thresholded";
if (!File.exists(tmp)) File.makeDirectory(tmp);

tmp = Input + "Channel2_Thresholded";
if (!File.exists(tmp)) File.makeDirectory(tmp);

tmp = Input + "Common_Area";
if (!File.exists(tmp)) File.makeDirectory(tmp);

tmp = Input + "Channel1_Measurements";
if (!File.exists(tmp)) File.makeDirectory(tmp);

tmp = Input + "Channel2_Measurements";
if (!File.exists(tmp)) File.makeDirectory(tmp);

Output1 = Input + "Channel1_Thresholded/"; ///// EDITABLE
print(Output1);
Output2 = Input + "Channel2_Thresholded/"; ///// EDITABLE
print(Output2);
Output3 = Input + "Common_Area/"; ///// EDITABLE
print(Output3);
Output4 = Input + "Channel1_Measurements/"; ///// EDITABLE
print(Output4);
Output5 = Input + "Channel2_Measurements/"; ///// EDITABLE
print(Output5);

if (isOpen("ROI Manager")) {
     	selectWindow("ROI Manager");
     	run("Close");
  		}

		clear=roiManager("count");
		if(clear==0)	{
		}	else	{
		cleararray=newArray(clear);
		for(l=0; l<clear;l++)	{
		cleararray[l]=l;
		}
		roiManager("Select", cleararray);
		roiManager("Delete");
		}

for(j=0;j<list.length;j++){ 
		open(list[j]);
		//run("Smooth");
		run("Sharpen");
		print(list[j]);
		title= getTitle();
		name = substring(title,0,lengthOf(title)-4);
		rename(name);
		selectWindow(name);
		roiManager("Add");
		//setBackgroundColor(0, 0, 0);
		//run("Clear Outside");
		run("Split Channels");
		
		selectWindow("C1-"+name);
		run("Enhance Contrast", "saturated=0.35 process_all");
		roiManager("Select",0);
		run("Measure");
		MEAN = getResult("Mean");
		STDEV = getResult("StdDev");
		CUTOFF1 = MEAN + STDEV + MOD1;
		if(CUTOFF1 > 255){CUTOFF1 = 1;} else{}
		run("Clear Results");
		print(CUTOFF1);
		
		selectWindow("C2-"+name);
		run("Enhance Contrast", "saturated=0.35 process_all");
		roiManager("Select",0);
		run("Measure");
		MEAN = getResult("Mean");
		STDEV = getResult("StdDev");
		CUTOFF2 = MEAN + STDEV + MOD2;
		if(CUTOFF2 > 255){CUTOFF2 = 1;} else{}
		run("Clear Results");
		print(CUTOFF2);

		if (isOpen("ROI Manager")) {
     	selectWindow("ROI Manager");
     	run("Close");
  		}

		clear=roiManager("count");
		if(clear==0)	{
		}	else	{
		cleararray=newArray(clear);
		for(l=0; l<clear;l++)	{
		cleararray[l]=l;
		}
		roiManager("Select", cleararray);
		roiManager("Delete");
		}
		
		// Start with Channel 1
		selectWindow("C1-"+name);
		setOption("BlackBackground", true);
		setThreshold(CUTOFF1, 255);
		run("Convert to Mask");
		run("Close-");
		run("Fill Holes");
		run("Watershed");
		//waitForUser;
		
		if (isOpen("ROI Manager")) {
     	selectWindow("ROI Manager");
     	run("Close");
  		}

		clear=roiManager("count");
		if(clear==0)	{
		}	else	{
			cleararray=newArray(clear);
		for(l=0; l<clear;l++)	{
			cleararray[l]=l;
		}
		roiManager("Select", cleararray);
		roiManager("Delete");
		}
		
		selectWindow("C1-"+name);
		run("Analyze Particles...", "size=2-8000 circularity=0.01-1.00 exclude add"); // May have to adjust this for Images derived from 60x confocal
		//run("Analyze Particles...", "size=0-20000 circularity=0.01-1.00 exclude add");
		selectWindow("C1-"+name);
		run("Measure");
		selectWindow("Results");
		Res_out1 = Output1 + name + ".txt";
		saveAs(Res_out1);
		run("Clear Results");
		selectWindow("C1-"+name);
		Res_out5 = Output1 + name + ".tif";
		saveAs(Res_out5);
		
		cell=roiManager("count");
		
		for(n=0;n<cell;n++){
			
			selectWindow("C1-"+name);
			roiManager("Select",n);
			run("Measure");
		}
		
		selectWindow("Results");
		Res_out4 = Output4 + name + ".txt";
		saveAs(Res_out4);
		run("Clear Results");
		
		
		if (isOpen("ROI Manager")) {
     	selectWindow("ROI Manager");
     	run("Close");
  		}

		clear=roiManager("count");
		if(clear==0)	{
		}	else	{
		cleararray=newArray(clear);
		for(l=0; l<clear;l++)	{
		cleararray[l]=l;
		}
		roiManager("Select", cleararray);
		roiManager("Delete");
		}
		
	
		// Now do Channel2	
		selectWindow("C2-"+name);
		setOption("BlackBackground", true);
		setThreshold(CUTOFF2, 255);
		run("Convert to Mask");
		run("Close-");
		run("Fill Holes");
		run("Watershed");
		//waitForUser;
		
		if (isOpen("ROI Manager")) {
     	selectWindow("ROI Manager");
     	run("Close");
  		}
		run("Analyze Particles...", "size=2-8000 circularity=0.01-1.00 exclude add");
		selectWindow("C2-"+name);
		run("Measure");
		selectWindow("Results");
		Res_out2 = Output2 + name + ".txt";
		saveAs(Res_out2);
		run("Clear Results");
		selectWindow("C2-"+name);
		Res_out6 = Output2 + name + ".tif";
		saveAs(Res_out6);
		
		cell=roiManager("count");
		
		for(n=0;n<cell;n++){
			
			selectWindow("C2-"+name);
			roiManager("Select",n);
			run("Measure");
		}
		
		selectWindow("Results");
		Res_out5 = Output5 + name + ".txt";
		saveAs(Res_out5);
		run("Clear Results");
		
		
		if (isOpen("ROI Manager")) {
     	selectWindow("ROI Manager");
     	run("Close");
  		}

		clear=roiManager("count");
		if(clear==0)	{
		}	else	{
		cleararray=newArray(clear);
		for(l=0; l<clear;l++)	{
		cleararray[l]=l;
		}
		roiManager("Select", cleararray);
		roiManager("Delete");
		}
		
		imageCalculator("AND create", "C1-"+name,"C2-"+name);
		
		selectWindow("Result of " + "C1-" + name);
		run("Analyze Particles...", "pixel exclude add");
		selectWindow("Result of " + "C1-" + name);
		rename(name);
		selectWindow(name);
		run("Convert to Mask");
		selectWindow(name);
		run("Measure");
		selectWindow("Results");
		Res_out3 = Output3 + name + ".txt";
		saveAs(Res_out3);
		
		if (isOpen("ROI Manager")) {
     	selectWindow("ROI Manager");
     	run("Close");
  		}

		clear=roiManager("count");
		if(clear==0)	{
		}	else	{
		cleararray=newArray(clear);
		for(l=0; l<clear;l++)	{
		cleararray[l]=l;
		}
		roiManager("Select", cleararray);
		roiManager("Delete");
		}
		
		run("Clear Results");
		run("Close All");
		
		if (isOpen("Log")) {
         selectWindow("Log");
         run("Close" );
    	}	
}
print("DONE");
