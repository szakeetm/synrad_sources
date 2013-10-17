/*
  File:        File.h
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

#ifndef FILERWH
#define FILERWH

#include <stdio.h>
#include "Types.h"
#include <string>

#define READ_BUFFSIZE 4096

// -------------------------------------------------------

class Error {

public:
  Error(char *message);
  const char *GetMsg();

private:
  char msg[1024];

};

// -------------------------------------------------------

class FileUtils {

public:
  // Utils functions
  static int Exist(char *fileName);
  static char *GetPath(char *fileName);
};

// -------------------------------------------------------

class FileReader {

public:
  // Constructor/Destructor
  FileReader(char *fileName);
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
  void Write(char *s);

private:

  FILE *file;
  char fileName[2048];

};
  
std::string SplitFilename (const std::string& str);
std::string SplitPath(const std::string& str);

#endif /* FILERWH */


