#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <Windows.h>
#include <io.h>
#include <direct.h>

std::string exec(std::string command);
std::string exec_str(const char* cmd);
BOOL Exist(std::string fileName);
BOOL Exist(const char *fileName);
std::string GetPath(const std::string& str);

int main(int argc,char* argv[]) {
	//__debugbreak();
	char key;
	std::string result;
	if (argc<3) {
		std::cout<<"Usage: compress.exe FILE_TO_COMPRESS RENAMED_FILENAME {include_file1 include_file2 ...}";
		ShowWindow( GetConsoleWindow(), SW_RESTORE );
		std::cin>>key;
		return 0;
	}
	std::string command;
	std::string fileName;
	std::string fileNameWith7z;
	std::string fileNameGeometry;
	/*memcpy(fileName,argv[1],strlen(argv[1])*sizeof(char));
	fileName[strlen(argv[1])*sizeof(char)]='\0';*/
	fileName = argv[1];
	std::cout<<"\nargv0: "<<argv[0];
	std::cout<<"\nargv1: "<<argv[1];
	std::cout<<"\nargv2: "<<argv[2];
	//sprintf_s(fileNameWith7z,"%s7z",fileName);
	fileNameWith7z = fileName + "7z";
	if (!Exist("7za.exe")) {
		printf("\n7za.exe not found. Cannot compress.\n");
			std::cin>>key;
			return 0;
	}
	char *dir;
	/*dir = strrchr(fileName,'\\');
	memcpy(fileNameGeometry,fileName,sizeof(char)*(dir-fileName));
	fileNameGeometry[dir-fileName]=NULL;
	sprintf_s(fileNameGeometry,"%s\\%s",fileNameGeometry,argv[2]);*/
	fileNameGeometry = GetPath(fileName) + "\\" + argv[2];
	//sprintf_s(command,"move \"%s\" \"%s\"",fileName,fileNameGeometry);
	command = "move \"" + fileName + "\" \"" + fileNameGeometry + "\"";
	result=exec(command);
	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );
	//delete destination file (empty archive)
	//sprintf_s(command,"del \"%s\"",fileNameWith7z);
	command = "del \""+fileNameWith7z+"\"";
	result = exec(command); 
	//sprintf_s(command,"cmd /C \"pushd \"%s\"&&7za.exe u -t7z \"%s\" \"%s\"",CWD,fileNameWith7z,fileNameGeometry);
	command = "cmd /C \"pushd \"";
	command+=CWD;
	command+="\"&&7za.exe u -t7z \"" + fileNameWith7z + "\" \"" + fileNameGeometry + "\"";
	for (int i=3;i<argc;i++) { //include files
		BOOL duplicate=false;
		for (int j=3;!duplicate && j<i;j++) { //check for duplicate include files
			if (strcmp(argv[i],argv[j])==0)
				duplicate=true;
		}
		if (!duplicate) {
			command += " \"";
			command += argv[i];
			command += "\""; //add as new input file
		}
	}
	command+="&&popd\"";
	std::cout<<"\nCommand: "<<command<<"\n\nStarting compression...\nYou can continue using Synrad while compressing.\n";
	result=exec(command);
	size_t found;
	//printf("\nresult: %s\n",result);
	found=result.find("Everything is Ok");
	if (found!=std::string::npos) {
		printf("\nCompression seems legit. Deleting GEO file.");
		remove(fileNameGeometry.c_str());
		return 0;
	}
	//printf("\nresult: %s\n",result);
	ShowWindow( GetConsoleWindow(), SW_RESTORE );
	//sprintf_s(command,"move \"%s\" \"%s\"",fileNameGeometry,fileName);
	command = "move \""+fileNameGeometry+"\" \""+fileName+"\"";
	result=exec(command);
	printf("\nSomething went wrong during the compression, read above. GEO file kept."
		"\nType any letter and press Enter to exit\n");
	std::cin>>key;
	return 0;
}

std::string exec(std::string command) {
	return exec_str(command.c_str());
}

std::string exec_str(const char* cmd) {
    FILE* pipe = _popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
        if(fgets(buffer, 128, pipe) != NULL)
                result += buffer;
		        printf(buffer);
    }
	result=result+'0';
    _pclose(pipe);
    return result;
}

BOOL Exist(std::string fileName) {
	return Exist(fileName.c_str());
}

BOOL Exist(const char *fileName) {

	if (FILE *file = fopen(fileName, "r")) {
		fclose(file);
		return TRUE;
	}
	return FALSE;
}

std::string GetPath(const std::string& str)
{
	size_t found = str.find_last_of("/\\");
	if (found == std::string::npos) return str; //not found
	else return str.substr(0, found);
}