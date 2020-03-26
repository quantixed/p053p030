/*
 * A macro to analyse knocksideways results in mitotic images
 * a pre and post knocksideways image with two channels (green and red)
 */


#@ File (label = "Input directory:", style = "directory") input

run("Clear Results");

if (File.exists(input+"/results.csv")){
	File.delete(input+"/results.csv");
}
f = File.open(input+"/results.csv");
print(f, "file,type,time,measure");
processFolder(input);

File.close(f);

function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	counter = 0;
	for (i = 0; i < list.length; i++) {
		if(endsWith(list[i], ".tif") && lengthOf(list[i]) > 15) {
			if(counter > 1) exit("This macro works on two images per folder only");
			if(counter == 0) pre = list[i];
			if(counter == 1) {
				post = list[i];
				processPair(input, pre, post);
			}
			counter += 1;
		}
	}
}

function processPair(input, pre, post) {
	open(input + File.separator + pre);
	preName = getTitle();
	open(input + File.separator + post);
	postName = getTitle();
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
	// centre the images - find the centres
	setTool("multipoint");
	selectWindow(preName);
	waitForUser("Use the multipoint tool to select the centre of the cell");
	getSelectionCoordinates(xpoints, ypoints);
	xPre = xpoints[xpoints.length - 1];
	yPre = ypoints[ypoints.length - 1];
	run("Select None");
	// post image
	setTool("multipoint");
	selectWindow(postName);
	waitForUser("Use the multipoint tool to select the centre of the cell");
	getSelectionCoordinates(xpoints, ypoints);
	xPost = xpoints[xpoints.length - 1];
	yPost = ypoints[ypoints.length - 1];
	run("Select None");
	// difference is in integers so can be used for translation
	xDiff = xPre - xPost;
	yDiff = yPre - yPost;
//	print(xPre + " " + yPre + " " + xPost + " " + yPost);
	run("Translate...", "x=" + xDiff + " y=" + yDiff + " interpolation=None stack");
	selectWindow(preName);
	run("Duplicate...", "title=pregreen duplicate channels=" + d2s(greenChannel,0));
	run("Grays");
	selectWindow(postName);
	run("Duplicate...", "title=postgreen duplicate channels=" + d2s(greenChannel,0));
	run("Grays");
	run("Merge Channels...", "c1=postgreen c2=pregreen create");
	mergeImage = "Composite";
	getROIs(mergeImage);
	processFile(preName,greenChannel,redChannel,0);
	processFile(postName,greenChannel,redChannel,1);
}

function getROIs(mergeImage){
	run("Set Measurements...", "area mean standard min integrated redirect=None decimal=3");
	// uncheck black background in options
	run("Options...", "iterations=1 count=1");
	// clear ROI Manager if there was some junk leftover
	roiManager("reset");
	
	selectWindow(mergeImage);
	getDimensions(xx, yy, cc, ss, ff);
	run("Select None");
	// first ROI is background
	setTool("rectangle");
	waitForUser("Use the rectangle select tool to select the background");
	roiManager("add");
	selectWindow(mergeImage);
	run("Select None");
	// now specify the spindle
	setTool("multipoint");
	waitForUser("Use the multipoint tool to select three points in the spindle");
	multipoint2rects(xx,yy);
	selectWindow(mergeImage);
	run("Select None");
	// now specify the cytoplasm
	setTool("multipoint");
	waitForUser("Use the multipoint tool to select three points in the cytoplasm");
	multipoint2rects(xx,yy);
	selectWindow(mergeImage);
	run("Select None");
	// now specify the mitochondria
	setTool("multipoint");
	waitForUser("Use the multipoint tool to select three points in the mitochondria");
	multipoint2rects(xx,yy);
	selectWindow(mergeImage);
	run("Select None");
	// now specify the cell
	setTool("freehand");
	waitForUser("Use the freehand select tool to draw around the cell");
	roiManager("Add");
	run("Select None");
	close();
}

function multipoint2rects(xx, yy)	{
	boxSize = 16;
	getSelectionCoordinates(xpoints, ypoints);
	newImage("TempImage", "8-bit black", xx, yy, 1);
	setColor(255, 255, 255);
	selectWindow("TempImage");
	for (i = 0; i < xpoints.length; i++) {
		fillRect(xpoints[i] - (boxSize / 2), ypoints[i] - (boxSize / 2), boxSize, boxSize);
	}
	run("Create Selection");
	run("Make Inverse");
	roiManager("Add");
	close();
}

function processFile(imgName, greenChannel, redChannel, counter) {
	res1 = getData(imgName, counter);

	print(f, imgName + "," + "spindle green," + counter + "," + res1[0]);
	print(f, imgName + "," + "cyto green," + counter + "," + res1[1]);
	print(f, imgName + "," + "mito green," + counter + "," + res1[2]);
	print(f, imgName + "," + "cell green," + counter + "," + res1[3]);
	print(f, imgName + "," + "spindle red," + counter + "," + res1[4]);
	print(f, imgName + "," + "cyto red," + counter + "," + res1[5]);
	print(f, imgName + "," + "mito red," + counter + "," + res1[6]);
	print(f, imgName + "," + "cell red," + counter + "," + res1[7]);
}


function getData(imgName, counter) {
	selectWindow(imgName);
	dir = getDirectory("image"); 
	results = newArray(8);
	for (i=0; i<5; i++)	{
		selectWindow(imgName);
		roiManager("Select", i);
		Stack.setPosition(greenChannel, 1, 1);
		run("Measure");
	}
	for (i=0; i<5; i++)	{
		selectWindow(imgName);
		roiManager("Select", i);
		Stack.setPosition(redChannel, 1, 1);
		run("Measure");
	}
	results[0] = getResult("Mean", 1) - getResult("Mean", 0);
	results[1] = getResult("Mean", 2) - getResult("Mean", 0);
	results[2] = getResult("Mean", 3) - getResult("Mean", 0);
	results[3] = getResult("Mean", 4) - getResult("Mean", 0);
	results[4] = getResult("Mean", 6) - getResult("Mean", 5);
	results[5] = getResult("Mean", 7) - getResult("Mean", 5);
	results[6] = getResult("Mean", 8) - getResult("Mean", 5);
	results[7] = getResult("Mean", 9) - getResult("Mean", 5);
	saveAs("results", dir + File.separator + imgName + "_measurements.csv");
	run("Clear Results");
	if(counter == 1) {
		roiManager("Save", dir + File.separator + "roiset.zip");
		run("Close All");
	}
	return results;
}