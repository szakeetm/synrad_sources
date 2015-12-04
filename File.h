/*
  File:        File.h
  Description: File management class
  Program:     SynRad
  Author:      R. KERSEVAN / M ADY
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

#ifndef FILERWH
#define FILERWH

#include <stdio.h>
#include "Types.h"
#include <string>

#define READ_BUFFSIZE 4096

// -------------------------------------------------------

class Error {

public:
  Error(const char *message);
  const char *GetMsg();

private:
  char msg[1024];

};

// -------------------------------------------------------

class FileUtils {

public:
  // Utils functions
	static int Exist(std::string fileName);
	static int Exist(const char *fileName);
	static std::string GetPath(const std::string &str); //Extracts string up to to last "\" (inlcuding "\"). If no path found, returns empty string
	static std::string GetFilename(const std::string &str); //Extracts string after the last "\"
	static std::string GetExtension(const std::string &str); //Extracts string after the last "."
};

// -------------------------------------------------------

class FileReader {

public:
  // Constructor/Destructor
	FileReader(std::string fileName) :FileReader(fileName.c_str()){};
	FileReader(const char *fileName);
	~FileReader();

  char *GetName();

  // Read function
  int IsEof();
  int IsEol();
  char *ReadLine();
  char *ReadString();
  llong ReadLLong();
  int ReadInt();
  double ReadDouble();
  void ReadKeyword(char *keyword);
  char *ReadWord();
  void JumpSection(char *end);
  void SeekStart();
  bool SeekFor(char *keyword);
  bool SeekForChar(char *c);
  bool wasLineEnd;

  Error MakeError(char *msg);
  int GetCurrentLine();

  void JumpComment();
  void JumpControlChars();
private:

  void RefillBuffer();
  char ReadChar();
  
  
  FILE *file;
  int curLine;
  char fileName[2048];
  char readBuffer[READ_BUFFSIZE];
  int  nbLeft;
  int  buffPos;
  int  isEof;
  char CurrentChar;
};

// -------------------------------------------------------

class FileWriter {

public:
  // Constructor/Destructor
  FileWriter(char *fileName);
  ~FileWriter();

  char *GetName();

  // Write function
  void WriteLLong(const llong &v,char *sep=NULL);
  void WriteInt(const int &v,char *sep=NULL);
  void WriteDouble(const double &v,char *sep=NULL);
  void Write(const char *s);
  void Write(std::string str);

private:

  FILE *file;
  char fileName[2048];

};

#endif /* FILERWH */


