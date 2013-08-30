/*
  File:        File.cpp
  Description: File management class
  Program:     SynRad
  Author:      R. KERSEVAN / M SZAKACS
  Copyright:   E.S.R.F / CERN

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "File.h"
#include <string.h>

#ifdef WIN32
#include <direct.h>
#include <io.h>
#include <iostream>
#include "GLApp\GLMessageBox.h"
#endif

#define MAX_WORD_LENGTH 65536 //expected length of the longest line

// -------------------------------------------------
// Error class
// -------------------------------------------------

Error::Error(char *message) {
  strncpy(msg,message,1023);
  msg[1023]=0;
}

// -------------------------------------------------

const char *Error::GetMsg() {
  return msg;
}

// -------------------------------------------------
// FileUtils class
// -------------------------------------------------

int FileUtils::Exist(char *fileName) {

#ifdef WIN32
  struct _finddata_t seqfile;
  intptr_t h;

  if( (h=_findfirst( fileName , &seqfile )) != -1L ) {
	  _findclose(h);
	}

  return (h != -1L);
#else
 // TODO
 return FALSE;
#endif
}

// -------------------------------------------------

char *FileUtils::GetPath(char *fileName) {

  static char tmp[512];
  strcpy(tmp,fileName);

  char  *p = strrchr(tmp,'\\');
  if(!p) p = strrchr(tmp,'/');
  if(p)  *p=0;
  return tmp;

}

// -------------------------------------------------
// FileReader class
// -------------------------------------------------

FileReader::FileReader(char *fileName) {

  file = fopen(fileName,"r");
  if(!file) throw Error("Cannot open file for reading");
  curLine = 1;
  strcpy(this->fileName,fileName);
  isEof = 0;
  CurrentChar = ' ';
  RefillBuffer();

}

// -------------------------------------------------

void FileReader::RefillBuffer() {

  if( !isEof ) {
    nbLeft = (int)fread(readBuffer,1,READ_BUFFSIZE,file);
    isEof = (nbLeft==0);
  }
  buffPos = 0;

}

// -------------------------------------------------

char FileReader::ReadChar() {

  if( buffPos>=nbLeft ) RefillBuffer();

  if( !isEof ) {
    CurrentChar = readBuffer[buffPos++];
	if(CurrentChar=='\n') {
		curLine++;
		wasLineEnd=true;
	}
  } else {
    CurrentChar = 0;
  }
  return CurrentChar;
}

// -------------------------------------------------

int FileReader::IsEof() {

  JumpControlChars();
  return isEof;

}

// -------------------------------------------------

int FileReader::GetCurrentLine() {

  return curLine;

}

// -------------------------------------------------

char *FileReader::GetName() {
  return fileName;
}

// -------------------------------------------------

FileReader::~FileReader() {
  fclose(file);
}

// -------------------------------------------------

Error FileReader::MakeError(char *msg) {
  static char ret[4096];
  sprintf(ret,"%s (line %d)",msg,curLine);
  return Error(ret);
}

// -------------------------------------------------

int FileReader::ReadInt() {

  int ret;
  char *w = ReadWord();
  if( sscanf(w,"%d",&ret)<=0 ) throw Error(MakeError("Wrong integer format"));
  return ret;

}

// -------------------------------------------------

llong FileReader::ReadLLong() {

  llong ret;
  char *w = ReadWord();
  if( sscanf(w,"%I64d",&ret)<=0 ) throw Error(MakeError("Wrong integer64 format"));
  return ret;

}

// -------------------------------------------------

void FileReader::SeekStart() {
  fseek(file, 0L, SEEK_SET);
  isEof = 0;
  curLine = 1;
  RefillBuffer();
  CurrentChar = ' ';
}

// -------------------------------------------------

void FileReader::JumpSection(char *end) {

  char *w = ReadWord();
  while(strcmp(w,end)!=0) w=ReadWord();

}

// -------------------------------------------------

void FileReader::ReadKeyword(char *keyword) {

  char *w = ReadWord();
  if( strcmp(w,keyword)!=0 ) {
	  char msg[200];
	  sprintf(msg,"Unexpected keyword in FileReader::ReadkeyWord()\n\"%s\" expected, \"%s\" found.",keyword,w);
	  throw Error(MakeError(msg));
  }
}

// -------------------------------------------------

double FileReader::ReadDouble() {

  double ret;
  char *w = ReadWord();
  if( sscanf(w,"%lf",&ret)<=0 ) {
	  throw Error(MakeError("Wrong double format"));
  }
  return ret;

}

// -------------------------------------------------
bool FileReader::SeekFor(char *keyword) {
	char *w;
	int i=0;
	bool notFound;
	do {
		w=ReadLine();
		notFound=((strcmp(w,keyword)!=0)&&(!isEof));
	} while (notFound);
	return (isEof)?false:true;
}

bool FileReader::SeekForChar(char *c) {
	char w;
	int i=0;
	bool notFound;
	do {
		w=this->ReadChar();
		notFound=(w!=(*c))&&(!isEof);
	} while (notFound);
	return (isEof)?false:true;
}

// -------------------------------------------------

char *FileReader::ReadLine() {

  static char retWord[MAX_WORD_LENGTH];
  int len = 0;
  
  /* Jump space and comments */
  JumpControlChars();

  
	  while( CurrentChar>=32 && len<MAX_WORD_LENGTH)
	  {
		retWord[len]=CurrentChar;
		len++;
		ReadChar();
	  }

	  if( len>=MAX_WORD_LENGTH )
	  {
		throw MakeError("Line too long");
		return NULL;
	  }


  retWord[len]='\0';
  return retWord;

}

// -------------------------------------------------

char *FileReader::ReadString()
{
  static char retStr[MAX_WORD_LENGTH];

  char *w = ReadWord();
  int len = (int)strlen(w);
  if( w[0]=='"' ) {
    strncpy(retStr,w+1,len-2);
    retStr[len-2]=0;
  } else {
    strcpy(retStr,w);
  }

  return retStr;

}

void FileReader::JumpComment() {
	JumpControlChars();
	if (CurrentChar=='{') {
		while (!isEof && CurrentChar!='}' )
			ReadChar();
		ReadChar();
	}
}

// -------------------------------------------------

void FileReader::JumpControlChars() {

  // Jump spaces and control characters
  while( !isEof && CurrentChar<=32 )
    ReadChar();

}

int FileReader::IsEol() {
	return CurrentChar=='\n';
}

// -------------------------------------------------

char *FileReader::ReadWord() {

  static char retWord[MAX_WORD_LENGTH];
  int len = 0;

  /* Jump space and comments */
  JumpControlChars();

  /* Treat special character */
  if( CurrentChar==':' || CurrentChar=='{' || CurrentChar=='}' || CurrentChar==',' )
  {
    retWord[0]=CurrentChar;
    retWord[1]='\0';
    ReadChar();
    return retWord;
  }

  /* Treat string */
  if( CurrentChar=='"' ) {
    retWord[len]=CurrentChar;
    len++;
    ReadChar();
    while( CurrentChar!='"' && CurrentChar!=NULL && CurrentChar!='\n' && len<MAX_WORD_LENGTH) {
      retWord[len]=CurrentChar;
      len++;
      ReadChar();
    }
    if( len>=MAX_WORD_LENGTH )
      throw MakeError("String too long");
    if( CurrentChar==NULL || CurrentChar=='\n')
      throw MakeError("String not ended");
    retWord[len]=CurrentChar;
    len++;
    ReadChar();
    retWord[len]='\0';
    return retWord;
  }

  /* Treat other word */
  while( CurrentChar>32 && CurrentChar!=':' && CurrentChar!='{'
    && CurrentChar!='}' && CurrentChar!=',' && len<MAX_WORD_LENGTH)
  {
    retWord[len]=CurrentChar;
    len++;
    ReadChar();
  }

  if( len>=MAX_WORD_LENGTH )
  {
    throw MakeError("String too long");
    return NULL;
  }

  retWord[len]='\0';
  return retWord;

}

// -------------------------------------------------
// FileWriter class
// -------------------------------------------------

FileWriter::FileWriter(char *fileName) {

  file = fopen(fileName,"w");
  if(!file) {
	  char tmp[256];
	  sprintf(tmp,"Cannot open file for writing %s",fileName);
	  throw Error(tmp);
  }
  strcpy(this->fileName,fileName);

}

// -------------------------------------------------

char *FileWriter::GetName() {
  return fileName;
}

// -------------------------------------------------

FileWriter::~FileWriter() {
  fclose(file);
}

void FileWriter::WriteInt(int v,char *sep) {
  if( !fprintf(file,"%d",v) )
    throw Error("Error while writing to file");
  if(sep) fprintf(file,"%s",sep);
}

void FileWriter::WriteLLong(llong v,char *sep) {
#ifdef WIN32
  if( !fprintf(file,"%I64d",v) )
    throw Error("Error while writing to file");
#else
  throw Error("FileWriter::WriteLLong() not implemented");
#endif
  if(sep) fprintf(file,"%s",sep);
}

void FileWriter::WriteDouble(double v,char *sep) {
  //if(v>=0.0) fprintf(file," ");
	fprintf(file," ");
  if( !fprintf(file,"%.14E",v) )
    throw Error("Error while writing to file");
  if(sep) fprintf(file,"%s",sep);
}

void FileWriter::Write(char *s) {
  if( !fprintf(file,"%s",s) )
    throw Error("Error while writing to file");
}

std::string SplitFilename (const std::string& str)
{
  size_t found;
  found=str.find_last_of("/\\");
  return str.substr(found+1);
}

std::string SplitPath (const std::string& str)
{
  size_t found;
  found=str.find_last_of("/\\");
  return str.substr(0,found);
}