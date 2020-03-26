/*
 * Original macro by Erick Martins Ratamero
 * Modified by Stephen Royle
 */


#@ File (label = "Input directory:", style = "directory") input

run("Clear Results");

if (File.exists(input+"/results.csv")){
	File.delete(input+"/results.csv");
}
f = File.open(input+"/results.csv");
print(f, "file,type,measure");
processFolder(input);

File.close(f);

function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	counter = 0;
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], ".tif") && lengthOf(list[i]) > 15) {
			if(counter > 1) exit("This macro works on two images per folder only");
			processFile(input, list[i], counter);
			counter += 1;
		}
	}
}

function processFile(input, file, counter) {
//	print("Processing: " + input + File.separator + file);
	open(input + File.separator + file);
	getROIs(file, counter);
	res1 = getData(file, counter);

	print(f, file + "," + "spindle green," + res1[0]);
	print(f, file + "," + "cell green," + res1[1]);
	print(f, file + "," + "spindle red," + res1[2]);
	print(f, file + "," + "cell red," + res1[3]);
	print(f, file + "," + "ratio green," + res1[4]);
	print(f, file + "," + "ratio red," + res1[5]);

	close(input + File.separator + file);
}


function getROIs(file, counter){
	title = getTitle();
	dir = getDirectory("image");
	run("Set Measurements...", "area mean standard min integrated redirect=None decimal=3");
	// clear ROI Manager if there was some junk leftover
	if(counter == 0) roiManager("reset");
	
	// first ROI is background
	setTool("rectangle");
	if(counter == 1) roiManager("Select", 0);
	waitForUser("Use the rectangle select tool to select the background");
	roiManager("Add");
	run("Select None");
	// now specify the cell
	setTool("freehand");
	if(counter == 1) roiManager("Select", 1);
	waitForUser("Use the freehand select tool to select the cell");
	roiManager("Add");
	run("Create Mask");
	selectWindow("Mask");
	rename("maskbig");
	// back to image
	selectWindow(file);
	run("Select None");
	// now specify two spindle halves
	setTool("freehand");
	if(counter == 1) roiManager("Select", 2);
	waitForUser("Use the freehand select tool to select the spindle. Use shift to select both half spindles");
	roiManager("Add");
	// make new ROI which encompasses both spindle halves and chromosomes
	run("Convex Hull");
	run("Create Mask");
	imageCalculator("XOR create", "maskbig","Mask");
	selectWindow("Result of maskbig");
	run("Create Selection");
	// this ROI is the cell minus the spindle + chromosomes
	roiManager("Add");
	selectWindow("Mask");
	close();
	selectWindow("maskbig");
	close();
	selectWindow("Result of maskbig");
	close();
}


function getData(file, counter){
	// which channel is red and which one green?
	// presume 1 = green and 2 = red
	greenChannel = 1;
	redChannel = 2;
	Stack.setPosition(1, 1, 1);
	getLut(reds, greens, blues);
	if(reds[255] > greens[255]) {
		greenChannel = 2;
		redChannel = 1;
	}
	
	dir = getDirectory("image"); 
	results = newArray(6);
	offset = counter * 4;

	roiManager("Select", 0 + offset);
	Stack.setPosition(greenChannel, 1, 1);
	run("Measure");
	roiManager("Select", 2 + offset);
	Stack.setPosition(greenChannel, 1, 1);
	run("Measure");
	roiManager("Select", 3 + offset);
	Stack.setPosition(greenChannel, 1, 1);
	run("Measure");

	roiManager("Select", 0 + offset);
	Stack.setPosition(redChannel, 1, 1);
	run("Measure");
	roiManager("Select", 2 + offset);
	Stack.setPosition(redChannel, 1, 1);
	run("Measure");
	roiManager("Select", 3 + offset);
	Stack.setPosition(redChannel, 1, 1);
	run("Measure");
	
	results[0] = getResult("Mean", 1) - getResult("Mean", 0);
	results[1] = getResult("Mean", 2) - getResult("Mean", 0);
	results[2] = getResult("Mean", 4) - getResult("Mean", 3);
	results[3] = getResult("Mean", 5) - getResult("Mean", 3);
	results[4] = results[0] / results[1];
	results[5] = results[2] / results[3];
	saveAs("results", dir + File.separator + file + "_measurements.csv");
	roiManager("Save", dir + File.separator + file + "_roiset.zip");
	if(counter == 1) roiManager("reset");
	run("Close All");
	run("Clear Results");
	return results;
}