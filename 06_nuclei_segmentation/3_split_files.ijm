// Final step of processing before using in CellProfiler
// This script will take .tif and stardist.tif files and split to Ch and Z single files
// Each image sequence will be saved in a seperate folder for each embryo

// Files can then be used in CellProfiler

setBatchMode(true);
print("\\Clear");
dir = getDirectory("Choose file dir");
print("dir: "+dir);
list = getFileList(dir);
print("Number o files: "+list.length);
fs = File.separator;
print(fs)

save_dir=dir+fs+"split_files";
File.makeDirectory(save_dir);
print("Saving files in: "+save_dir);

for(i=0;i<list.length;i++) {
	if(endsWith(list[i], ".tif") == true && endsWith(list[i], "_stardist.tif") == false) {
		print("File nr: "+i+1);
		run("Bio-Formats Importer", "open=["+dir+fs+list[i]+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
		title = File.nameWithoutExtension();
		print("Image title: "+title);
		run("Bio-Formats Exporter", "save=["+save_dir+fs+title+fs+title+".tif] write_each_z_section write_each_channel compression=Uncompressed");
		run("Close All");
		print("done");
	};
	if(endsWith(list[i], "_stardist.tif") == true) {
		print("File nr: "+i+1);
		print("Stardist File");
		run("Bio-Formats Importer", "open=["+dir+fs+list[i]+"] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
		title = File.nameWithoutExtension();
		titleshort = replace(title, "_stardist", "");
		print("Image title: "+title);
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		run("Bio-Formats Exporter", "save=["+save_dir+fs+titleshort+fs+title+".tif] write_each_z_section compression=Uncompressed");
		run("Close All");
		print("done");
	};
	showProgress(-i/list.length);
};
print("All done");
setBatchMode(false);