// This macro will open all series of a lif file
// and save as .tif multichannel z stack
// Use to process with Stardist and CellProfiler scripts
// If you uncomment the last line, it will also convert the stack to images, saving all C and Z in its own folder


path = File.openDialog("Select a File");
output = getDirectory("Select directory");
fs = File.separator;

run("Bio-Formats Macro Extensions");
Ext.setId(path);
Ext.getSeriesCount(seriesCount);

for (s=1; s<=seriesCount; s++) {
	run("Bio-Formats Importer", "open=&path autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_"+s);

    title = getTitle;
    title = replace(title, ".lif", "");
    title = replace(title, " ", ""); 
    title = replace(title, "/", "_"); //
    print(title);
	saveAs("Tiff", output+title+".tif");

	close();
}

//setBatchMode(true); 

//for (i=0;i<nImages;i++) {
//        selectImage(i+1); 
//        title = getTitle;
//        title = replace(title, ".lif", "");
//        title = replace(title, " ", ""); 
        //title = replace(title, ".", "");
//        print(title);
        //save_dir = output+fs+title;
		//File.makeDirectory(save_dir);
		//print("Saving files in: "+save_dir);
        //name = getInfo("image.filename");
		//print(name);
		//run("Bio-Formats Exporter", "save=["+save_dir+fs+title+".tif] write_each_z_section write_each_channel compression=Uncompressed");
//		saveAs("Tiff", output+title+".tif");
//	};


//etBatchMode(false); 
//run("Close All");