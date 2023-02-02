dir = getDirectory("Please select a directory");
print(dir)
filenames = getFileList(dir);

WrongPattern1 = ".*Z[0-9]{1,2}\\.tif$";
WrongPattern2 = ".*Z[0-9]{1,2}_C.*\\.tif$";

for(i = 0; i < filenames.length; i++) {
	newFileName = filenames[i];
	while (matches(newFileName, WrongPattern1) | matches(newFileName, WrongPattern2)) {
		print("Wrong Pattern Matched");
		newFileName = replace(newFileName, "Z", "Z0");
	} 
	File.rename(dir+filenames[i], dir+newFileName);
}

print("Done")
