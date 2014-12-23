#ifndef _CONFIGFILE_H_
#define _CONFIGFILE_H_

/*

* Copyright (c) 2007, Carnegie Mellon University
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of Carnegie Mellon University, nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  Code author: Jernej Barbic
  CMU, 2005-2007
  Version 1.0

  A class to parse text configuration files. 
  See configFile-example.cpp and configFile-example.config for examples. 
  It is by far easier to look at the example first, than to proceed
  straight to reading the class interface.

  Supported types for option entries:
  int (integer)
  bool (boolean, true/false)
  float (single precision floating point)
  double (double precision floating point)
  char* (a C-style string)

  Configuration files are text files. The format for each option entry is:
  *<option name>
  <option value>

  Example:
  *temperature
  36.6

  Blank lines are ignored. 
  A line is a comment if it begins with a "#".
  The options can appear in arbirary order. They need not follow the order
  in which they were created via addOption/addOptionOptional .

*/

#include <vector>
#include <string>

class ConfigFile
{
public:
 
  ConfigFile();
  ~ConfigFile();

  // === routines to define the valid option entries for your configuration file ===
  // each option entry is either mandatory, or optional (in which case you need to provide a default value)
  // option names are case IN-sensitive (i.e., upper/lower case does not matter)
  // "destLocation" will be overwritten with the value read from the configuration file when "parseOptions" is called (or default value will be used in case of optional option entries if a particular configuration file did not provide the entry)

  // routines to specify mandatory option entries
  // if "parseOptions" does not find a mandatory option in a particular configuration file, it will exit with a non-zero code
  int addOption(const char * optionName, int * destLocation);
  int addOption(const char * optionName, bool * destLocation);
  int addOption(const char * optionName, float * destLocation);
  int addOption(const char * optionName, double * destLocation);
  int addOption(const char * optionName, char * destLocation); // for strings, you must pass a pointer to a pre-allocated string buffer (and not a pointer to a (char*) pointer)

  // routines to specify option entries that are optional
  // if not specified in a particular config file, the value will default to the given default value
  template<class T>
  int addOptionOptional(const char * optionName, T * destLocation, T defaultValue);
  int addOptionOptional(const char * optionName, char * destLocation, const char * defaultValue);

  // === after you have specified your option entries with addOption and/or addOptionOptional, call "parseOptions" to open a particular configuration file and load the option values to their destination locations.

  int parseOptions(const char * filename, int verbose=1); // returns 0 on success, and a non-zero value on failure
  int parseOptions(FILE * fin, int verbose=1); // makes it possible to parse from a file stream

  // after calling "parseOptions", you can print out the values of all options, to see the values that were read from the configuration file
  void printOptions();

  // (optional) set a stopping string (when it is encountered in a file, parsing will stop); default: **EOF
  // if the stopping string does not appear in a file, it will be parsed to the end
  void setStoppingString(const char * stoppingString);

  // you can disable printing out warnings (default: enabled)
  void suppressWarnings(int suppressWarnings_) { this->suppressWarnings_ = suppressWarnings_; }

protected:
  std::vector<std::string> optionNames;
  std::vector<int> optionTypes;
  std::vector<void*> destLocations;

  std::vector<bool> optionSet;
  
  int seekOption(const char * optionName); // returns -1 if not found

  template<class T>
  int addOptionHelper(const char * optionName, T * destLocation);

  int getTypeSize(int type);
  void getTypeFormatSpecifier(int type, char * fsp);

  void upperCase(char * s); // converts a string to upper case
  void removeTrailingCharacters(char * s, char ch); // removes (one or more) characters 'ch' from the end of the string

  char stoppingString[32];

  int suppressWarnings_;
};

#endif

