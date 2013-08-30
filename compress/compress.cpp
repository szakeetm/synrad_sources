#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <Windows.h>
#include <io.h>
#include <direct.h>

std::string exec(char* cmd);
int FileExists(char *fileName);

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
	char command[20000];
	char fileName[2048];
	char fileNameWith7z[2048];
	char fileNameGeometry[2048];
	memcpy(fileName,argv[1],strlen(argv[1])*sizeof(char));
	fileName[strlen(argv[1])*sizeof(char)]='\0';
	std::cout<<"\nargv0: "<<argv[0];
	std::cout<<"\nargv1: "<<argv[1];
	std::cout<<"\nargv2: "<<argv[2];
	sprintf_s(fileNameWith7z,"%s7z",fileName);
	if (!FileExists("7za.exe")) {
		printf("\n7za.exe not found. Cannot compress.\n");
			std::cin>>key;
			return 0;
	}
	char *dir;
	dir = strrchr(fileName,'\\');
	memcpy(fileNameGeometry,fileName,sizeof(char)*(dir-fileName));
	fileNameGeometry[dir-fileName]=NULL;
	sprintf_s(fileNameGeometry,"%s\\%s",fileNameGeometry,argv[2]);
	sprintf_s(command,"move \"%s\" \"%s\"",fileName,fileNameGeometry);
	result=exec(command);
	char CWD [MAX_PATH];
	_getcwd( CWD, MAX_PATH );
	//delete destination file (empty archive)
	sprintf_s(command,"del \"%s\"",fileNameWith7z);
	result=exec(command);
	sprintf_s(command,"cmd /C \"pushd \"%s\"&&7za.exe u -t7z \"%s\" \"%s\"",CWD,fileNameWith7z,fileNameGeometry);
	for (int i=3;i<argc;i++) { //include files
		BOOL duplicate=false;
		for (int j=3;!duplicate && j<i;j++) { //check for duplicate include files
			if (strcmp(argv[i],argv[j])==0)
				duplicate=true;
		}
		if (!duplicate) sprintf_s(command,"%s \"%s\"",command,argv[i]); //add as new input file
	}
	sprintf_s(command,"%s&&popd\"",command);
	std::cout<<"\nCommand: "<<command<<"\n\nStarting compression...\nYou can continue using Synrad while compressing.\n";
	result=exec(command);
	size_t found;
	//printf("\nresult: %s\n",result);
	found=result.find("Everything is Ok");
	if (found!=std::string::npos) {
		printf("\nCompression seems legit. Deleting GEO file.");
		remove(fileNameGeometry);
		return 0;
	}
	//printf("\nresult: %s\n",result);
	ShowWindow( GetConsoleWindow(), SW_RESTORE );
	sprintf_s(command,"move \"%s\" \"%s\"",fileNameGeometry,fileName);
	result=exec(command);
	printf("\nSomething went wrong during the compression, read above. GEO file kept."
		"\nType any letter and press Enter to exit\n");
	std::cin>>key;
	return 0;
}

std::string exec(char* cmd) {
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

int FileExists(char *fileName) {

  struct _finddata_t seqfile;
  intptr_t h;

  if( (h=_findfirst( fileName , &seqfile )) != -1L ) {
	  _findclose(h);
	}

  return (h != -1L);
}