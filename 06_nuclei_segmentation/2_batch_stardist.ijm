// This script will run Stardist segmentation on .tif multichannel Z stack images

/*
 * Macro template to process multiple images in a folder
 */
 // Run LIF series exporter FIRST
 // Input files are .tif multichannel stack
 // Will save Stardist in the same file
 // After this run Split_stack

#@ File (label = "Input directory", style = "directory") input // Input of .tif files
#@ File (label = "Output directory", style = "directory") output // Output should be the same folder for subsequent scripts to work
#@ String (label = "File suffix", value = ".tif") suffix

// See also Process_Folder.py for a version of this code
// in the Python scripting language.
print("\\Clear");
processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
// Note : make sure not additional subfolders and files!
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	print("Processing: " + input + File.separator + file);
	open(input+File.separator + file);
	filename_short = File.nameWithoutExtension;
	save_file = filename_short + "_stardist.tif";
	save_path = output + File.separator + save_file;
	run("Reduce Dimensionality...", "slices"); // Single Z slice = DAPI CH1, would need to alter if segmenting on anything else
	// title = getTitle;
	// run("Split Channels");
	// selectWindow("C1-"+title);
	run("Re-order Hyperstack ..." , "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]"); // to trick stardist doing pseudo 3D must invert Z and t
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'"+file+"', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', 'probThresh':'0.5', 'nmsThresh':'0.5', 'outputType':'Both', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	selectWindow("Label Image");
	print("Saving to in reality: " + save_path);
	save(save_path);
	print("Saving to: " + output);
	run("Close All");
}
