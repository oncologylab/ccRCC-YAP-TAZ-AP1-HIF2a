// Input is RAW .Tiff Files in Micron Scale
// Goal of this script is to isolate nuclei and to output individual nuclei, without background using DAPI as a marker, so puncta information can later be obtained from other channels.
run("Close All");
run("Set Measurements...", "area mean min standard display redirect=None decimal=3");

Directories = "/Volumes/.../Directories.txt"; // plain text file, tab separated, containing name(s) of folder(s) used for each experiment in which IF (2 channel + DAPI) was performed, In our case, we used the files directly from a portable harddrive
d = File.openAsString(Directories);
d = split(d, "\r\n");
Array.print(d);

Groups = "/Volumes/.../Filenames.txt"; // plain text file, tab separated, containing prefix names of conditions used to name image files. The repetitions of these should reflect
g = File.openAsString(Groups);
g = split(g, "\r\n");
Array.print(g);

for(f=0;f<d.length;f++){
	
	print("Start: " + d[f]); 
	print("Group: " + g[f]);
	Input = d[f];
	list = getFileList(Input);
	Array.print(list);
	
	MainDir = "/Users/.../Final_Confocal_Images/"+ g[f] +"/"; ///// EDITABLE
	if (!File.exists(MainDir)) File.makeDirectory(MainDir);
	
	tmp = MainDir + "Combined_Images1";
	if (!File.exists(tmp)) File.makeDirectory(tmp);
	
	tmp1 = MainDir + "Nuclear_Area_Stats1";
	if (!File.exists(tmp1)) File.makeDirectory(tmp1);
		
	tmp2 = MainDir + "Ch01_Stats";
	if (!File.exists(tmp2)) File.makeDirectory(tmp2);
		
	tmp3 = MainDir + "Ch02_Stats";
	if (!File.exists(tmp3)) File.makeDirectory(tmp3);
	
	tmp4 = MainDir + "Nuclear_Images1";
	if (!File.exists(tmp4)) File.makeDirectory(tmp4);
		
	Output1 = "/Users/.../Final_Confocal_Images/" + g[f] + "/Combined_Images1/"; ///// EDITABLE
	Output2 = "/Users/.../Final_Confocal_Images/" + g[f] + "/Nuclear_Area_Stats1/"; ///// EDITABLE
	Output3 = "/Users/.../Final_Confocal_Images/" + g[f] + "/Ch01_Stats/"; ///// EDITABLE
	Output4 = "/Users/.../Final_Confocal_Images/" + g[f] + "/Ch02_Stats/"; ///// EDITABLE
	Output5 = "/Users/.../Final_Confocal_Images/" + g[f] + "/Nuclear_Images1/"; ///// EDITABLE
	
	
	b = File.openAsString(Input + "/" + "names.txt");
	b = split(b, "\r\n");
	
	for(s=0; s<b.length;s++){
		
		Samp_list = Array.concat(Samp_list, b[s]);
		sample = b[s];
		print(sample);
	
		
		files2=newArray(); 
		
		for(i=0;i<list.length;i++){ 
			if(endsWith(list[i],"ch00.tif")&&indexOf(list[i],sample)>=1){files2=Array.concat(files2,list[i]);} 
		}
		Array.sort(files2);
		Array.print(files2);
		Znum = lengthOf(files2);
		
		for(j=0;j<files2.length;j++){ 
			open(Input + "/" + files2[j]);
			print(files2[j]);
			run("8-bit");
			run("Enhance Contrast", "saturated=0.35");
			run("Subtract Background...", "rolling=200");
		}
		run("Images to Stack", "use");
		sliceMax = nSlices();
		print(sliceMax);
		if(sliceMax > 9){nchar=13;} else{nchar=12;}
		run("Z Project...", "projection=[Max Intensity]");
		
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
			
			//title= getTitle();
			//name = substring(title,0,lengthOf(title)-nchar);
			run("8-bit");
			run("Enhance Contrast", "saturated=0.35");
			setAutoThreshold("Default dark");
			setThreshold(25, 255); //////////////// EDITABLE ////////////// If nuclei are too faded, lower
			run("Convert to Mask");
			run("Fill Holes");
			run("Despeckle");
			run("Analyze Particles...", "size=90000-750000 pixel exclude add");
			// Add in if statement regarding 4 
			cell=roiManager("count");
			run("Close All");   
		    
		    for(j=0;j<files2.length;j++){ 
			open(Input + "/" + files2[j]);
			run("8-bit");
			//run("Enhance Contrast", "saturated=0.35");
			//run("Subtract Background...", "rolling=200");
			title= getTitle();
			print(title);
			}
		run("Images to Stack", "use");
		selectWindow("Stack");
		rename("DAPI");
		
		name = substring(title,0,lengthOf(title)-nchar);
		print(name);
		
		//Read in Channel1
		files3=newArray();
		
		for(i=0;i<list.length;i++){ 
			if(endsWith(list[i],"ch01.tif")&&indexOf(list[i],sample)>=1){files3=Array.concat(files3,list[i]);} 
			}
			
		for(j=0;j<files3.length;j++){ 
			open(Input + "/" + files3[j]);
			print(files3[j]);
			run("8-bit");
		}
		run("Images to Stack", "use");	 
		selectWindow("Stack");
		rename("Channel1");
		
		//Read in Channel2
		files4=newArray();
		
		for(i=0;i<list.length;i++){ 
			if(endsWith(list[i],"ch02.tif")&&indexOf(list[i],sample)>=1){files4=Array.concat(files4,list[i]);} 
			}
			
		for(j=0;j<files4.length;j++){ 
			open(Input + "/" + files4[j]);
			print(files4[j]);
			run("8-bit");
		}
		run("Images to Stack", "use");	 
		selectWindow("Stack");
		rename("Channel2");
		
		selectWindow("DAPI");
		run("Z Project...", "projection=[Max Intensity]");
		selectWindow("Channel1");
		run("Z Project...", "projection=[Max Intensity]");
		selectWindow("Channel2");
		run("Z Project...", "projection=[Max Intensity]");
		selectWindow("DAPI");
		close();
		selectWindow("Channel1");
		close();
		selectWindow("Channel2");
		close();
		
		// Merge Channels for image to use in publication
		run("Merge Channels...", "c1=[MAX_Channel1] c2=[MAX_Channel2] create keep");
		
		
		
			for(n=0;n<cell;n++){
			print("Cell:"+n);
			selectWindow("MAX_DAPI");
			roiManager("Select",n);
			run("Duplicate...", "duplicate");
			rename(name+"_"+n);
			run("Measure");
			run("Clear Outside");
			tmp = getTitle();
			selectWindow(tmp);
			Res_out5 = Output5 + name + "_" + n + ".tif";
			saveAs(Res_out5);
			close();
			selectWindow("Results");
			Res_out1 = Output2 + name + "_" + n + ".txt";
			saveAs(Res_out1);
			run("Clear Results");
			
			//open up the other two channels
			//Ch01
			selectWindow("MAX_Channel1");
			roiManager("Select",n);
			run("Duplicate...", "duplicate");
			rename(tmp);
			run("Measure");
			close();
			selectWindow("Results");
			Res_out3 = Output3 + name + "_" + n + ".txt";
			saveAs(Res_out3);
			run("Clear Results");
			//selectWindow(tmp);
			//rename("Channel1");
			
			
			// Ch02
			selectWindow("MAX_Channel2");
			roiManager("Select",n);
			run("Duplicate...", "duplicate");
			rename(tmp);
			run("Measure");
			close();
			selectWindow("Results");
			Res_out4 = Output4 + name + "_" + n + ".txt";
			saveAs(Res_out4);
			run("Clear Results");
			//selectWindow(tmp);
			//rename("Channel2");
			
			
			// Save merged image to use in publication
			selectWindow("Composite");
			roiManager("Select",n);
			run("Duplicate...", "duplicate");
			rename(tmp);
			run("Clear Outside");
			run("Enhance Contrast", "saturated=0.15 process_all");
			selectWindow(tmp);
			Res_out2 = Output1 + name + "_" + n + ".tif";
			saveAs(Res_out2);
			close();
			
			if (isOpen("Log")) {
         selectWindow("Log");
         run("Close" );
    			}	
			}
			
			run("Close All");
		
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
			
			 
			
			
	}
print("Sample Complete");	
	
	}
print("Done");