/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2015 Robert Baron, Nandita Mukhopadhyay,
  Lee Almasy, Mark Schroeder, William P. Mulvihill,
  Daniel E. Weeks, and University of Pittsburgh

  This file is part of the Mega2 program, which is free software; you
  can redistribute it and/or modify it under the terms of the GNU
  General Public License as published by the Free Software Foundation;
  either version 3 of the License, or (at your option) any later
  version.

  Mega2 is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

  For further information contact:
      Daniel E. Weeks
      e-mail: weeks@pitt.edu

===========================================================================
*/

#ifndef BATCH_H
#define BATCH_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cctype>
#include <algorithm>
#include <functional>
#include <vector>
#include <map>
#include <errno.h>

#include "common.h"
#include "analysis.h"
#include "analysis_ext.h"
#include "batch_input.h"
#include "reorder_loci_ext.h"


using namespace std;


//
// If this is defined file and directory checks will be disabled...
#if 0
#define DISABLE_FILE_CHECKS
#endif

//
// If these are defined a BatchFileException will not be thrown in the associated circumstance...
//#define DISABLE_EXCEPTION_IF_KEYWORD_DEPRECATED

#define DISABLE_EXCEPTION_ANALYSIS_DEPENDS_ON_SUBANALYSIS


static inline string toUpper(string str)
{
    for (size_t pos = 0; str[pos] != '\0'; ++pos)
        str[pos] = (char)toupper(str[pos]);
    return str;
}

/**
 Batch file processing error handling.
 
 The 'BatchFileException' is thrown if an error is found while processing a batch file.
 The exception has associated with one of these enumerations as well as a more descriptive string.
 */
typedef enum  {
  BATCHFILE_UNKNOWN_EXCEPTION, BATCHFILE_BATCHITEM_LINE_INVALID, BATCHFILE_BATCHITEM_DEPENDENCY_NOT_MET, BATCHFILE_BATCHITEM_REQUIRED, BATCHFILE_BATCHITEM_INVALID_KEYWORD, BATCHFILE_BATCHITEM_NOT_FOUND, BATCHFILE_BATCHITEM_SUBANALYSIS_MANDITORY, BATCHFILE_DUPLICATE_ENTRY, BATCHFILE_FILE_HAS_BEEN_PROCESSED, BATCHFILE_FILE_NOT_FOUND, BATCHFILE_VALUE_NOT_INTEGER, BATCHFILE_VALUE_NOT_INTEGER_GT0, BATCHFILE_VALUE_NOT_INTEGER_GE0, BATCHFILE_VALUE_NOT_INTEGER_VECTOR, BATCHFILE_VALUE_NOT_INTEGER_VECTOR_ETERM, BATCHFILE_VALUE_NOT_DOUBLE, BATCHFILE_VALUE_NOT_DOUBLE_GE0, BATCHFILE_VALUE_NOT_DOUBLE_VECTOR, BATCHFILE_VALUE_NOT_BOOLEAN, BATCHFILE_DEPRECATED_ENTRY, BATCHFILE_INTERNAL_ERROR, BATCHFILE_INVALID_METHOD, BATCHFILE_ANALYSIS_OPTION_UNKNOWN, BATCHFILE_ANALYSIS_OPTION_DISABLED, BATCHFILE_ANALYSIS_OPTION_CHANGE, BATCHFILE_NO_VALUE_FOR_ENTRY, BATCHFILE_VALUE_NOT_IN_PERMITTED_LIST, BATCHFILE_BATCHITEM_PLINK_BAD_PARAMETER, BATCHFILE_BATCHITEM_INVALID_DIRECTORY, BATCHFILE_BATCHITEM_FILE_DOES_NOT_EXIST
} BatchFileExceptionType;

/**
 Thrown when some fatal error occurs in batch file processing...
 */
class BatchFileException {
private:
    string s;
    BatchFileExceptionType type;
    
public:
    BatchFileException(std::string ss) : s(ss) { type = BATCHFILE_UNKNOWN_EXCEPTION; }
    BatchFileException(std::string ss, BatchFileExceptionType type) : s(ss), type(type) {}
    ~BatchFileException() {}
    
    const char* what() const { return s.c_str(); }
    BatchFileExceptionType getType() const { return type; }
};


class BatchFile; // defined here beacuse it is a friend of BatchItem

/**
   @brief BatchItem

   This is the base class that contains information about a batch item.
   The base class holds book keeping information that is common to all
   batch items. Derived classes provide more specialized functionality
   and informaiton.
 
   The old style batch file had key value pairs (e.g., key = value).
   This implementation alows for the notion of a section for grouping
   keys in much the same way that the old style Windows .INI files did.
   This class only records the string associated with the value, and some
   book keeping informaiton.
 
   Every possible BatchItem is instantiated by BatchFile (a friend).
   They are protected so that no one else can instantiate them.
*/
class BatchItem {
 private:
  // Taken together the keyword/section provides a unique identifier for the
  // batch item. For backwards compatibility the empty section (e.g., "") is used.
  string keyword;
  string section;
    
  // BIG NOTE: The 'valueString' should be pulled from the base class since it will not be valid after the value is reset.
  // The raw string after being trimed.
  string valueString;
    
  // book keeping information about the batch item.
    
  // Some items are always required, others based on what other items have been specified.
  // Because of this no provision is made for setting this in the creator.
  bool isRequired;
    
  // Whether the item was read from the batch file...
  bool wasReadFromBatchFile;
    
  // Whe the value has been changed from the one read from the batch file,
  // or it has been assigned and never read...
  bool valueWasReplaced;
    
  // This is typically assigned on instantiation...
  bool hasDefaultValue;
    
  bool isDeprecated;
    
  // Must have at least one virtual method to the class to be considered polymorphic.
  // Without this the dynamic casts do not work...
  virtual void unused() {};
    
 protected:
  // There is not enough information in the base class to instantiate one in and
  // of itself, so the creator is protected...
 BatchItem(string keyword, string section) : keyword(keyword), section(section) {
    valueString = "";
    isRequired = false;
    wasReadFromBatchFile = false;
    valueWasReplaced = false;
    hasDefaultValue = false;
    isDeprecated = false;
    // for case insensitive compares...
    this->keyword = keyword;
    this->section = section;
  }
    
  // setters... Only used by methods in the derived classes...
  void setWasReadFromBatchFile() { wasReadFromBatchFile = true; }
  void resetWasReadFromBatchFile() { wasReadFromBatchFile = false; }
  void setHasDefaultValue() { hasDefaultValue = true; }
  void setIsDeprecated() { isDeprecated = true; }
  void setValueWasReplaced() { valueWasReplaced = true; }
    
  // The method 'setValue()' should be used ONLY when reading a batch file,
  // that is why it is protected. If you want to modify the value that has
  // been read from the batch file, or assign a value when a batch file
  // was never present, you should use the 'replaceValue()' method which
  // is defined for each different subclass of BatchItem. The 'replace'
  // method is not virtual since the parameter that it takes will be different
  // for each subclass of BatchItem.
    
  // BIG NOTE: The 'valueString' should be pulled from the base class since it will
  // not be valid after the value is reset.
  virtual void setValue(string valueString_) { this->valueString = valueString_; }
    
 public:
  virtual ~BatchItem() {}

  // setters...

  void setIsRequired() { isRequired = true; }
        
  // getters...
        
  string& getKeyword() { return keyword; }
  char const *getKeywordCstr() const { return keyword.c_str(); }
        
  string& getSection() { return section; }
        
  string& getValueString() { return valueString; }
  char const *getValueStringCstr() const { return valueString.c_str(); }
        
  bool itemRead() const { return wasReadFromBatchFile; }
  const int itemReadIntValue() const { return (itemRead() ? 1 : 0); }
  bool getIsRequired() const { return isRequired; }
  bool getHasDefaultValue() const { return hasDefaultValue; }
  bool getIsDeprecated() const { return isDeprecated; }
  bool getValueWasReplaced() const { return valueWasReplaced; }
        
    
  // This won't work for batchItems because it is a list of pointers...
  bool operator<(const BatchItem& rhs) const {
    return (const string)section < (const string)rhs.section && (const string)keyword < (const string)rhs.keyword;
  }
        
  // This allows us to do "cout << instanceBatchItem;" and have the appropriate thing written...
  // This is pure virtual, so EACH derived class must implement it...
  friend ofstream& operator<< (ofstream& os, const BatchItem& bi);


  // This could be an abstract method, but it's not clear what to do if you are
  // asking for the value string for a deprecated item?!
  virtual const char *getValueCstr() const { return ""; }

  // This allows us to do "cout << instanceBatchItem;" and have the appropriate thing written...
  virtual void print(ofstream& os) const {
    os << getValueCstr();
  }

  // Determine if there is a dependency on this BatchItem that is not being
  // met because the argument has not been read. If so throw an exception.
  void checkDependencyOn(BatchItem bi) {
    if (itemRead() && !bi.itemRead()) {
      stringstream ss;
      ss << "Batch Item '" << getKeyword();
      if (getSection() != "") ss << "' (Section '" << getSection() << "')";
      ss << "' is dependent on '" << bi.getKeyword();
      if (bi.getSection() != "") ss << "' (Section '" << bi.getSection() << "')";
      ss << "' which is not specified.",
	throw BatchFileException(ss.str(), BATCHFILE_BATCHITEM_DEPENDENCY_NOT_MET);
    }
  }
        
  // friends...
        
  // So the predicate can see our private variables...
  //    friend class BatchItemPointerCompare;
  friend class BatchItemKeywordEqual;
  friend class BatchItemSectionEqual;
  friend class BatchItemSectionKeywordEqual;
        
  // So that it can see use 'wasReadFromBatchFile()'...
  friend class BatchFile;
};
        
inline ofstream& operator<<(ofstream& os, const BatchItem& bi){
  bi.print(os);
  return os;
}
#if 0
struct BatchItemPointerCompare {
  bool operator() (BatchItem *lhs, BatchItem *rhs) {
    return lhs->getSection() < rhs->getSection() && lhs->getKeyword() < rhs->getKeyword();
  }
};
#endif
        
/**
   Thes classes define an operator (Boolean) which is used to determine if the BatchItem
   is associated with the default section ""/keyword string. This is done with case insensitive compares.
*/
class BatchItemKeywordEqual : public binary_function <BatchItem *, string, bool> {
 public:
  bool operator() (BatchItem *lhs, const string keyword) const {
    return lhs->getSection() == "" && toUpper(lhs->getKeyword()) == toUpper(keyword);
  }
};
class BatchItemSectionEqual : public binary_function <BatchItem *, string, bool> {
 public:
  bool operator() (BatchItem *lhs, const string section) const {
    return lhs->getSection() == section;
  }
};
/**
   Thes classes define an operator (Boolean) which is used to determine if the BatchItem
   is associated with the section/keyword string which is packaged into a vector.
   This is done with case insensitive compares.
*/
class BatchItemSectionKeywordEqual : public binary_function <BatchItem *, vector<string>, bool> {
 public:
  bool operator() (BatchItem *lhs, const vector<string> sk) const {
    return toUpper(lhs->getSection()) == toUpper(sk[0]) && toUpper(lhs->getKeyword()) == toUpper(sk[1]);
  }
};
        
        
/**
   Here we create derived classes of the batch item. They follow the patern of the
   original batch input code including a: boolean (YORN), int, double (float), string,
   integer list, integer list with legal values, string list, string list with legal values,
   float list.
*/

class BatchItemInt : public BatchItem {
 protected:
  int value;
  int defaultValue;
            
 BatchItemInt(string keyword, string section) : BatchItem(keyword, section) { value = 0; }
 BatchItemInt(string keyword, string section, int defaultValue) : BatchItem(keyword, section) {
    setHasDefaultValue();
    value = this->defaultValue = defaultValue;
  }
            
  // see the documentation for BatchItem::setValue() for a description...
  void setValue(string valueString_) {
    // BIG NOTE: The 'valueString' should be pulled from the base class since it will not be valid after the value is reset.
    BatchItem::setValue(valueString_); // record the string in the base class.
    stringstream ss(valueString_);
    if ((ss >> value).fail()) {
      throw BatchFileException("Batch Item '" + getKeyword() + "' value '" + valueString_ + "' should be an Integer.",
			       BATCHFILE_VALUE_NOT_INTEGER);
    }
  }
            
 public:
  // see the documentation for BatchItem::setValue() for a description...
  void replaceValue(string valueString_) {
    setValue(valueString_);
    setValueWasReplaced();
  }
  void replaceValue(int value_) {
    this->value = value_;
    setValueWasReplaced();
  }
  int getValue() { return value; }
  virtual const char* getValueCstr() const {
    ostringstream os;
    os << value;
    return os.str().c_str();
  }
  int getDefaultValue() { return defaultValue; }
            
  // So that it can instantiate this item...
  friend class BatchFile;
};
        
        
class BatchItemIntOf : public BatchItemInt {
private:
    vector<int> valuesPermitted;
    
    void valueCheck() {
        if (find(valuesPermitted.begin(), valuesPermitted.end(), BatchItemInt::getValue()) == valuesPermitted.end()) {
            stringstream ss;
            ss << "Batch Item '" << getKeyword() << ", the value '" << BatchItemInt::getValue() << "' is not in the list of permitted values {";
            for (vector<int>::iterator it=valuesPermitted.begin(); it != valuesPermitted.end(); ++it) {
                ss << *it;
                if (it+1 != valuesPermitted.end()) ss << " ";
            }
            ss << "}";
            throw BatchFileException(ss.str(), BATCHFILE_VALUE_NOT_IN_PERMITTED_LIST);
        }
    }
    
protected:
    BatchItemIntOf(string keyword, string section, vector<int>valuesPermitted) : BatchItemInt(keyword, section) {
        this->valuesPermitted = valuesPermitted;
    }
    BatchItemIntOf(string keyword, string section, vector<int>valuesPermitted, int defaultValue) : BatchItemInt(keyword, section, defaultValue) {
        this->valuesPermitted = valuesPermitted;
    }
    
    // see the documentation for BatchItem::setValue() for a description...
    void setValue(string valueString_) {
        // BIG NOTE: The 'valueString_' should be pulled from the base class since it will not be valid after the value is reset.
        BatchItemInt::setValue(valueString_);
        // Now check to see if the value in the valuesPermitted vector is legal...
        valueCheck();
    }
    
public:
    // see the documentation for BatchItem::setValue() for a description...
    void replaceValue(int value_) {
        BatchItemInt::replaceValue(value_);
        valueCheck();
    }
    
    // So that it can instantiate this item...
    friend class BatchFile;
};
        
        
/**
   Same as BatchItemInt but with an extra check to make sure that the value is always > 0.
*/
class BatchItemIntGt0 : public BatchItemInt {
 protected:
 BatchItemIntGt0(string keyword, string section) : BatchItemInt(keyword, section) { }
 BatchItemIntGt0(string keyword, string section, int defaultValue) : BatchItemInt(keyword, section, defaultValue) { }
            
  void checkValue() {
    if (getValue() <= 0) {
      stringstream ss;
      ss << "Batch Item '" << getKeyword() << "' value '" << getValue() << "' should be an Integer greater than zero.";
      throw BatchFileException(ss.str(), BATCHFILE_VALUE_NOT_INTEGER_GT0);
    }
  }
            
  // see the documentation for BatchItem::setValue() for a description...
  void setValue(string valueString_) {
    BatchItemInt::setValue(valueString_);
    checkValue();
  }
            
 public:
  // see the documentation for BatchItem::setValue() for a description...
  void replaceValue(int value_) {
    BatchItemInt::replaceValue(value_);
    checkValue();
  }
            
  // so that it can instantiate this item...
  friend class BatchFile;
};

class BatchItemIntGe0 : public BatchItemInt {
 protected:
 BatchItemIntGe0(string keyword, string section) : BatchItemInt(keyword, section) { value = 1; }
 BatchItemIntGe0(string keyword, string section, int defaultValue) : BatchItemInt(keyword, section, defaultValue) { }
            
  void checkValue() {
    if (getValue() < 0) {
      stringstream ss;
      ss << "Batch Item '" << getKeyword() << "' value '" << getValue() << "' should be a positive Integer greater than or equal to zero.";
      throw BatchFileException(ss.str(), BATCHFILE_VALUE_NOT_INTEGER_GE0);
    }
  }
            
  // see the documentation for BatchItem::setValue() for a description...
  void setValue(string valueString_) {
    BatchItemInt::setValue(valueString_);
    checkValue();
  }
            
 public:
  // see the documentation for BatchItem::setValue() for a description...
  void replaceValue(int value_) {
    BatchItemInt::replaceValue(value_);
    checkValue();
  }
            
  // So that it can instantiate this item...
  friend class BatchFile;
};
        
class BatchItemIntVector : public BatchItem {
 private:
            
 protected:
  vector<int> value;
  vector<int> defaultValue;

 BatchItemIntVector(string keyword, string section) : BatchItem(keyword, section) { }
 BatchItemIntVector(string keyword, string section, vector<int> defaultValue) : BatchItem(keyword, section) {
    setHasDefaultValue();
    value = this->defaultValue = defaultValue;
  }
            
  // see the documentation for BatchItem::setValue() for a description...
  void setValue(string valueString_) {
    BatchItem::setValue(valueString_);
    stringstream ss(valueString_);
                
    // parse the integers and put them into the vector
    value.clear();
    for (int i; ss >> i;) { value.push_back(i); if (ss.eof()) break; }
                
    if (ss.fail()) {
      throw BatchFileException("Batch Item '" + getKeyword() + "' value '" + valueString_ + "' should be a Vector of Integers.",
			       BATCHFILE_VALUE_NOT_INTEGER_VECTOR);
    }
  }
            
 public:
  // see the documentation for BatchItem::setValue() for a description...
  void resetValue() { value = defaultValue; }
  void replaceValue(string valueString_) {
    setValue(valueString_);
    setValueWasReplaced();
  }
  void replaceValue(vector<int> value_) {
    this->value = value_;
    setValueWasReplaced();
  }
  // Process an array of values 'n' items in length...
  void replaceValue(int *values, int n) {
    value.clear();
    for (int i = 0; i < n; i++) value.push_back(values[i]);
    setValueWasReplaced();
  }
  void clear() {
    value.clear();
  }
  void appendToValue(int value_) {
    this->value.push_back(value_);
    setValueWasReplaced();
  }
  vector<int> getValue() { return value; }
  virtual const char *getValueCstr() const {
    ostringstream os;
    for (vector<int>::const_iterator it = value.begin(); it != value.end(); ++it) {
      os << *it;
      if (it+1 != value.end()) os << " ";
     }
    return (const char *)os.str().c_str();
  }
  const int getElement(const int i) { return value[(size_t)i]; }
  const size_t size() const { return value.size(); }
            
  // For backward compatability for C programs.
  // The space returned must be explicitly freed by the caller.
  // Alternatively, we could make one copy and free it in the distructor.
  int *getValueClist() const {
    int *clist = (int *)calloc(value.size(), sizeof(int));
    for (size_t i=0; i < value.size(); i++) clist[i] = value[i];
    return clist;
  }
            
  vector<int> getDefaultValue() { return defaultValue; }
            
  // So that it can instantiate this item...
  friend class BatchFile;
};

class BatchItemIntVectorChromosome : public BatchItemIntVector {
    
protected:
    BatchItemIntVectorChromosome(string keyword, string section) : BatchItemIntVector(keyword, section) { }
    BatchItemIntVectorChromosome(string keyword, string section, vector<int> defaultValue) : BatchItemIntVector(keyword, section, defaultValue) { }
    
    // see the documentation for BatchItem::setValue() for a description...
    void setValue(string valueString_) {
        BatchItem::setValue(valueString_);
        stringstream ss(valueString_);
        
        // parse the integers and put them into the vector
        for (string ch; ss >> ch;) {
            // STR_CHR will return -1 on error
            int i = STR_CHR(ch.c_str());
            if (i != -1) {
                value.push_back(i);
            } else {
                stringstream ss2;
                ss2 << "A chromosome must be a positive integer less than '" << lastautosome+4 << "' or one of the following: X, Y, XY, MT, U.\n";
                //ss << "A chromosome must be a positive integer or one of the following: X, Y, XY, MT, U.\n";
                throw BatchFileException(ss2.str(), BATCHFILE_VALUE_NOT_IN_PERMITTED_LIST);
            }
            if (ss.eof()) break;
        }
    }
    
    // So that it can instantiate this item...
    friend class BatchFile;
};

//
// There is some "cheating" done here. We only deal with the 'e' on the end if we are printing the value out or
// are taking it in from a string. Otherwise, a vector of ints is just fine.
class BatchItemIntVectorEterm : public BatchItemIntVector {
            
 protected:

 BatchItemIntVectorEterm(string keyword, string section) : BatchItemIntVector(keyword, section) { }
 BatchItemIntVectorEterm(string keyword, string section, vector<int> defaultValue) : BatchItemIntVector(keyword, section, defaultValue) { }
            
  // see the documentation for BatchItem::setValue() for a description...
  void setValue(string valueString_) {
    BatchItem::setValue(valueString_);
    stringstream ss(valueString_);
                
    // This should be a string of integers terminated by the 'e' character...
    for (int i; ss >> i;) {this->value.push_back(i); if (ss.eof()) break;}
                
    // If there wasn't an 'e' at the end we should have an else clause that complains
    // there no 'e' was found at the end... But, in this case, we will take a list of integers
    // without an 'e' at the end as well...
    if (ss.fail()) {
      ss.clear(); // the failure status
      ss.unget(); // the non-integer character that caused the failure
      string x;
      ss >> x; // get the remainder of the line as a string which and better be an "e"
      if (x != "e")
	throw BatchFileException("Batch Item '" + getKeyword() + "' value '" + valueString_ + "' should be a Vector of Integers terminated with the 'e' character.",
				 BATCHFILE_VALUE_NOT_INTEGER_VECTOR_ETERM);
    }        // parse the integers and put them into the vector
  }
            
 public:
  virtual const char *getValueCstr() const {
    ostringstream os;
    for (vector<int>::const_iterator it = value.begin(); it != value.end(); ++it) {
      os << *it;
      if (it+1 != value.end()) os << " ";
     }
    os << " e";
    return (const char *)os.str().c_str();
  }
            
  // So that it can instantiate this item...
  friend class BatchFile;
};
        
        
class BatchItemDouble : public BatchItem {
 private:
  double value;
  double defaultValue;
            
 protected:
 BatchItemDouble(string keyword, string section) : BatchItem(keyword, section) { value = 0.0; }
 BatchItemDouble(string keyword, string section, double defaultValue) : BatchItem(keyword, section) {
    setHasDefaultValue();
    value = this->defaultValue = defaultValue;
  }
            
  // see the documentation for BatchItem::setValue() for a description...
  void setValue(string valueString_) {
    BatchItem::setValue(valueString_);
    stringstream ss(valueString_);
    if ((ss >> value).fail()) {
      throw BatchFileException("Batch Item '" + getKeyword() + "' value '" + valueString_ + "' should be a Double.",
			       BATCHFILE_VALUE_NOT_DOUBLE);
    }
  }
  void setValue(double value_) {
    this->value = value_;
  }
            
 public:
  // see the documentation for BatchItem::setValue() for a description...
  void replaceValue(string valueString_) {
    setValue(valueString_);
    setValueWasReplaced();
  }
  void replaceValue(double value_) {
    this->value = value_;
    setValueWasReplaced();
  }
  double getValue() const { return value; }
  virtual const char* getValueCstr() const {
    ostringstream os;
    os << value;
    return os.str().c_str();
  }
  double getDefaultValue() const { return defaultValue; }
            
  // So that it can instantiate this item...
  friend class BatchFile;
};
        
        
class BatchItemDoubleGe0 : public BatchItemDouble {
 protected:
   BatchItemDoubleGe0(string keyword, string section) : BatchItemDouble(keyword, section) { }
   BatchItemDoubleGe0(string keyword, string section, double defaultValue) : BatchItemDouble(keyword, section, defaultValue) { }
            
  void checkValue() {
    if (getValue() < 0) {
      stringstream ss;
      ss << "Batch Item '" << getKeyword() << "' value '" << getValue() << "' should be an float greater than or equal to zero.";;
      throw BatchFileException(ss.str(), BATCHFILE_VALUE_NOT_DOUBLE_GE0);
    }
  }
            
  // see the documentation for BatchItem::setValue() for a description...
  void setValue(string valueString_) {
    BatchItemDouble::setValue(valueString_);
    checkValue();
  }
            
 public:
  // see the documentation for BatchItem::setValue() for a description...
  void replaceValue(int value_) {
    BatchItemDouble::replaceValue(value_);
    checkValue();
    setValueWasReplaced();
  }
            
  // So that it can instantiate this item...
  friend class BatchFile;
};
        
        
class BatchItemMissingQuant : public BatchItemDouble {
 protected:
 BatchItemMissingQuant(string keyword, string section) : BatchItemDouble(keyword, section) {}
 BatchItemMissingQuant(string keyword, string section, double defaultValue) : BatchItemDouble(keyword, section, defaultValue) {}
            
  // see the documentation for BatchItem::setValue() for a description...
  void setValue(string valueString_) {
    if (valueString_ == "NA") BatchItemDouble::replaceValue(QMISSING);
    else BatchItemDouble::setValue(valueString_);
  }
  void setValue(double value_) {
    BatchItemDouble::setValue(value_);
  }
            
 public:

  // So that it can instantiate this item...
  friend class BatchFile;
};
        
class BatchItemDoubleVector : public BatchItem {
 private:
  vector<double> value;
  vector<double> defaultValue;
            
 protected:
 BatchItemDoubleVector(string keyword, string section) : BatchItem(keyword, section) { }
 BatchItemDoubleVector(string keyword, string section, vector<double> defaultValue) : BatchItem(keyword, section) {
    setHasDefaultValue();
    value = this->defaultValue = defaultValue;
  }
            
  // see the documentation for BatchItem::setValue() for a description...
  void setValue(string valueString_) {
    BatchItem::setValue(valueString_);
    stringstream ss(valueString_);
                
    // parse the integers and put them into the vector
    for (double i; ss >> i; ) { value.push_back(i); if (ss.eof()) break; }
                
    if (ss.fail()) {
      throw BatchFileException("Batch Item '" + getKeyword() + "' value '" + valueString_ + "' should be a Vector of Doubles.",
			       BATCHFILE_VALUE_NOT_DOUBLE_VECTOR);
    }
  }
            
 public:
  // see the documentation for BatchItem::setValue() for a description...
  void replaceValue(string valueString_) {
    value.clear();
    setValue(valueString_);
    setValueWasReplaced();
  }
  void replaceValue(vector<double> value_) {
    this->value = value_;
    setValueWasReplaced();
  }
  // Process an array of values 'n' items in length...
  void replaceValue(double *values, int n) {
    value.clear();
    for (int i = 0; i < n; i++) value.push_back(values[i]);
    setValueWasReplaced();
  }
  void clear() {
    value.clear();
  }
  void appendToValue(double value_) {
    this->value.push_back(value_);
    setValueWasReplaced();
  }
  virtual const char *getValueCstr() const {
    ostringstream os;
    for (vector<double>::const_iterator it = value.begin(); it != value.end(); ++it) {
      os << *it;
      if (it+1 != value.end()) os << " ";
    }
    return (const char *)os.str().c_str();
  }            
  vector<double> getValue() { return value; }
  const double getElement(const int i) { return value[(size_t)i]; }
  const size_t size() const { return value.size(); }
            
  // For backward compatability for C programs.
  // The space returned must be explicitly freed by the caller.
  // Alternatively, we could make one copy and free it in the distructor.
  double *getValueClist() const {
    double *clist = (double *)calloc(value.size(), sizeof(double));
    for (size_t i=0; i < value.size(); i++) clist[i] = value[i];
    return clist;
  }

  vector<double> getDefaultValue() { return defaultValue; }
            
  // So that it can instantiate this item...
  friend class BatchFile;
};
        
        
class BatchItemBoolean : public BatchItem {
 private:
  bool value;
  bool defaultValue;
  static char mytoupper(char c) {  return (char)toupper(c); }
            
 protected:
 BatchItemBoolean(string keyword, string section) : BatchItem(keyword, section) { value = false; }
 BatchItemBoolean(string keyword, string section, bool defaultValue) : BatchItem(keyword, section) {
    setHasDefaultValue();
    value = this->defaultValue = defaultValue;
  }
            
  // see the documentation for BatchItem::setValue() for a description...
  void setValue(string valueString_) {
    BatchItem::setValue(valueString_);
    valueString_ = toUpper(valueString_);
    if (valueString_ == "YES" || valueString_ == "Y" || valueString_ == "TRUE" || valueString_ == "T" || valueString_ == "1") {
      value = true;
    } else if (valueString_ == "NO" || valueString_ == "N" || valueString_ == "FALSE" || valueString_ == "F" || valueString_ == "0") {
      value = false;
    } else {
      throw BatchFileException("Batch Item '" + getKeyword() + "' value '" + valueString_ + "' should be a Boolean (e.g.. yes, y, true, t, 1, no, n, false, f, 0).",
			       BATCHFILE_VALUE_NOT_DOUBLE);
    }
  }
            
 public:
  // see the documentation for BatchItem::setValue() for a description...
  void replaceValue(string valueString_) {
    setValue(valueString_);
    setValueWasReplaced();
  }
  void replaceValue(bool value_) {
    this->value = value_;
    setValueWasReplaced();
  }
  void replaceValueInt(int value_) {
    if (value_ == 1) this->value = true;
    else if (value_ == 0) this->value = false;
    else {
        stringstream ss;
        ss << value_;
        throw BatchFileException("Batch Item '" + getKeyword() + "' value '" + ss.str() + "' Integer boolean values must be either 1 or 0.",
                                 BATCHFILE_VALUE_NOT_BOOLEAN);
    }
    setValueWasReplaced();
  }
  void replaceValueYorN(int value_) {
    if (value_ == 'y' || value_ == 'Y') this->value = true;
    else if (value_ == 'n' || value_ == 'N') this->value = false;
    else {
	throw BatchFileException("Batch Item '" + getKeyword() + "' value must be either {y,Y,n,N}.",
                             BATCHFILE_VALUE_NOT_BOOLEAN);
    }
  }
  virtual const char *getValueCstr() const {
    return value ? "true" : "false";
  }

  bool getValue() { return value; }
            
  // For the "C" programs which use 1 as true and 0 as false...
  const int getIntValue() const { return value ? 1 : 0; }
  const char getCharValue() const { return value ? 'y' : 'n'; }
            
  bool getDefaultValue() { return defaultValue; }
            
  // So that it can instantiate this item...
  friend class BatchFile;
};
        
        
class BatchItemString : public BatchItem {
 private:
  string value;
  string defaultValue;
            
 protected:
 BatchItemString(string keyword, string section) : BatchItem(keyword, section) {}
 BatchItemString(string keyword, string section, string defaultValue) : BatchItem(keyword, section) {
    setHasDefaultValue();
    value = this->defaultValue = defaultValue;
  }
            
  // see the documentation for BatchItem::setValue() for a description...
  void setValue(string valueString_) {
    value = this->value = valueString_;
  }
            
 public:
  // see the documentation for BatchItem::setValue() for a description...
  void replaceValue(string valueString_) {
    value = valueString_;
    setValueWasReplaced();
  }
  string getValue() { return value; }
  virtual const char* getValueCstr() const { return value.c_str(); }

  string getDefaultValue() { return defaultValue; }
  const char* getDefaultValueCstr() { return defaultValue.c_str(); }

  const bool is_numeric() const {
    const char *value_ = this->value.c_str();
    char *end;
    // The converted value_ is thrown away, because we just want to know if it is valid...
    IgnoreValue(strtod((const char *)value_, &end));
    // See user_input.cpp:set_missing_quant_input() for an explaination of this test...
    return (strlen(value_) != 0 && strlen(end) == 0 && errno != ERANGE);
  }
            
  // So that it can instantiate this item...
  friend class BatchFile;
};
        

/**
   Class holds a string that can only be of the set 'valuesPermitted'...
*/
class BatchItemStringOf : public BatchItemString {
 private:
  vector<string> valuesPermitted;
            
  void valueCheck() {
    if (find(valuesPermitted.begin(), valuesPermitted.end(), BatchItemString::getValue()) == valuesPermitted.end()) {
      stringstream ss;
      ss << "Batch Item '" << getKeyword() << ", the value '" << BatchItemString::getValue() << "' is not in the list of permitted values {";
      for (vector<string>::iterator it=valuesPermitted.begin(); it != valuesPermitted.end(); ++it) {
	ss << *it;
	if (it+1 != valuesPermitted.end()) ss << " ";
      }
      ss << "}";
      throw BatchFileException(ss.str(), BATCHFILE_VALUE_NOT_IN_PERMITTED_LIST);
    }
  }
            
 protected:
 BatchItemStringOf(string keyword, string section, vector<string> valuesPermitted) : BatchItemString(keyword, section) {
    this->valuesPermitted = valuesPermitted;
  }
 BatchItemStringOf(string keyword, string section, vector<string> valuesPermitted, string defaultValue) : BatchItemString(keyword, section, defaultValue) {
    this->valuesPermitted = valuesPermitted;
  }
            
  // see the documentation for BatchItem::setValue() for a description...
  void setValue(string valueString_) {
    BatchItemString::setValue(valueString_);
    // Now check to see if the value in the valuesPermitted vector is legal...
    valueCheck();
  }
            
 public:
  void replaceValue(string value_) {
    BatchItemString::setValue(value_);
    // Now check to see if the value in the valuesPermitted vector is legal...
    valueCheck();
    setValueWasReplaced();
  }
            
  // So that it can instantiate this item...
  friend class BatchFile;
};

        
class BatchItemStringVector : public BatchItem {
 private:
  vector<string> defaultValue;
            
  // break the string into tokens separated by white space and
  // store them into the value vector...
  void tokenizeStringAndSetValue(string valueString_, string delimiters) {
    size_t startpos = 0;
    size_t pos = valueString_.find_first_of(delimiters, startpos);
    string token;
                
    while (string::npos != pos || string::npos != startpos) {
      token = valueString_.substr(startpos, pos - startpos);
      value.push_back(token.substr(0, token.length()));
                    
      startpos = valueString_.find_first_not_of(delimiters, pos);
      pos = valueString_.find_first_of(delimiters, startpos);
    }
                
    setWasReadFromBatchFile();
  }
            
 protected:
  vector<string> value;
            
  BatchItemStringVector(string keyword, string section) : BatchItem(keyword, section) { }
  BatchItemStringVector(string keyword, string section, vector<string> defaultValue) : BatchItem(keyword, section) {
    setHasDefaultValue();
    value = this->defaultValue = defaultValue;
  }
            
  void setValue(string valueString_) {
    BatchItem::setValue(valueString_);
    tokenizeStringAndSetValue(valueString_, " \t");
  }
            
 public:
  // see the documentation for BatchItem::setValue() for a description...
  void replaceValue(string valueString_) {
    value.clear();
    setValue(valueString_);
    setValueWasReplaced();
  }
  void replaceValue(vector<string> value_) {
    this->value = value_;
    setValueWasReplaced();
  }
  // Process an array of values 'n' items in length...
  void replaceValue(char *values[], int n) {
    value.clear();
    for (int i = 0; i < n; i++) value.push_back(values[i]);
    setValueWasReplaced();
  }
  virtual const char *getValueCstr() const {
    ostringstream os;
    for (vector<string>::const_iterator it = value.begin(); it != value.end(); ++it) {
      os << *it;
      if (it+1 != value.end()) os << " ";
    }
    return (const char *)os.str().c_str();
  }            
  vector<string> getValue() { return value; }
  const string getElement(const int i) { return value[(size_t)i]; }
  const char* getElementCstr(const int i) const { return value[(size_t)i].c_str(); }
  const size_t size() const { return value.size(); }
            
  // For backward compatability for C programs.
  // The space returned must be explicitly freed by the caller.
  // Alternatively, we could make one copy and free it in the distructor.
  const char **getValueClist() const {
    const char **clist = (const char **)calloc(value.size(), sizeof(char *));
    for (size_t i=0; i < value.size(); i++) clist[i] = this->getElementCstr((int)i);
    return clist;
  }

  vector<string> getDefaultValue() { return defaultValue; }
            
  // So that it can instantiate this item...
  friend class BatchFile;
};
        
//
// There is not a lot of value here at the moment... A method to open the file could be added...
// It mostly exists because there is a 'BatchItemFileInput'.
class BatchItemFile : public BatchItemString {
 protected:
 BatchItemFile(string keyword, string section) : BatchItemString(keyword, section) { }
 BatchItemFile(string keyword, string section, string defaultValue) : BatchItemString(keyword, section, defaultValue) { }

 public:
            
  // So that it can instantiate this item...
  friend class BatchFile;
};


class BatchItemFileInput : public BatchItemFile {
 protected:
 BatchItemFileInput(string keyword, string section) : BatchItemFile(keyword, section) { }
 BatchItemFileInput(string keyword, string section, string defaultValue) : BatchItemFile(keyword, section, defaultValue) { }

  void valueCheck() {
#ifndef DISABLE_FILE_CHECKS
    if (access(getValue().c_str(), F_OK) != 0) {
      throw BatchFileException("Batch Item '" + getKeyword() + "' the associated file '" + getValue() + "' is not found.",
			       BATCHFILE_BATCHITEM_FILE_DOES_NOT_EXIST);
    }
#endif /* DISABLE_FILE_CHECKS */
  }

  void setValue(string valueString_) {
    BatchItemFile::setValue(valueString_);
    valueCheck();
  }

 public:
  void replaceValue(string valueString_) {
    setValue(valueString_);
    setValueWasReplaced();
  }
            
  // So that it can instantiate this item...
  friend class BatchFile;
};


class BatchItemDirectory : public BatchItemString {
 protected:
 BatchItemDirectory(string keyword, string section) : BatchItemString(keyword, section) { }
 BatchItemDirectory(string keyword, string section, string defaultValue) : BatchItemString(keyword, section, defaultValue) { }

  void valueCheck(); // in batch_input.cpp

  void setValue(string valueString_) {
    BatchItemString::setValue(valueString_);
    valueCheck();
  }

 public:
  void replaceValue(string valueString_) {
    setValue(valueString_);
    setValueWasReplaced();
  }
  void replaceValueNoCheck(string valueString_) {
    BatchItemString::setValue(valueString_);
    setValueWasReplaced();
  }
            
  // So that it can instantiate this item...
  friend class BatchFile;
};
        
/**
   This class holds a pointer to the analysis type where it gets it's information about
   the particulars of the analysis.
 */        
class BatchItemAnalysis : public BatchItem {
 protected:
    static CLASS_ANALYSIS *value;

    BatchItemAnalysis(string keyword, string section) : BatchItem(keyword, section) { }
            
    void setValue(CLASS_ANALYSIS *value_) {
        this->value = value_;
    }
    void setValue(string valueString_) {
        BatchItem::setValue(valueString_);

        bool found = false;
	// This is a list of the "top level" analysis. That is to say, those analysis modes that
	// are offered to the user in the ANALYSIS MENU created in user_input.cpp::analysis_menu1()
        extern analysis_types analysis_list[];
        extern int count_analysis_list;

        // Here we check to see it the 'Analysis_Option' that we found in the batch file is assoiated
        // with a valid class by comparng the value to the CLASS_ANALYSIS::_name;
        for (int i = 0; i < count_analysis_list; i++) {
            analysis_types at = analysis_list[i];
            CLASS_ANALYSIS *a = at.analysis;
            const char *analysisName = a->_name;
            if (strcasecmp(analysisName, valueString_.c_str()) == 0) {
                found = true;
                value = a;
                break;
            }
        }
        if (!found) {
            throw BatchFileException("Batch Item '" + getKeyword() + "' the analysis option '" + valueString_ + "' is unknown.",
                                     BATCHFILE_ANALYSIS_OPTION_UNKNOWN);
        }
        setWasReadFromBatchFile();
    }
            
            
 public:
  // It does not make sense to replace the value of the analysis...
  //void replaceValue(vector<string> value) { }
  virtual const char *getValueCstr() const {
    return value->_name;
  }
            
  void replaceValue(string valueString_) {
    setValue(valueString_);
    setValueWasReplaced();
  }
  void replaceValue(CLASS_ANALYSIS *value_) {
    this->value = value_;
    setValueWasReplaced();
  }
            
  CLASS_ANALYSIS *getValue() { return value; }
            
  // So that it can instantiate this item...
  friend class BatchFile;
};


CLASS_ANALYSIS *BatchItemAnalysis::value = NULL;

        
class BatchItemDeprecated : public BatchItem {
 private:
  void  deprecatedException() {
    throw BatchFileException("Batch Item '" + getKeyword() + "' has been deprecated.",
			     BATCHFILE_DEPRECATED_ENTRY);
  }
            
 protected:
 BatchItemDeprecated(string keyword, string section) : BatchItem(keyword, section) {
    setIsDeprecated();
  }
            
 public:
  //string& getValueString() { deprecatedException(); }
  //string& getKeyword() { deprecatedException(); }
  //bool itemRead() { deprecatedException(); }
  //bool getHasDefaultValue() { deprecatedException(); }
            
  // So that it can instantiate this item...
  friend class BatchFile;
};
        
        
/**
   A batch file contains keyword value pairs one per line. A comment is permitted after
   the '#' character following a keyword value pair, or on a line without a keyword value pair.
   For example....
         
   # The comment character does not need to begin in the first column
   foo = 17 # this is the value of foo
         
   There is a notion of sections (named and generic) as there is in an old fassioned Windows
   .INI file. The sections are generally based on analysis types where the keywords appearing
   in those sections are specific to that analysis type or section name.
   For example...
         
   [ PLINK ]
   cm = yes # we will be reading genetic distances in centimorgans
   []
   Input_Locus_File = mylocus.txt
         
   Describes a "PLINK" section in which the keyword "cm" is defined, followed by the beginning
   of a "generic" section in which the keyword "Input_Locus_File" is defined. The batch file
   begins in the generic section for backward compatibility.
         
   This class contains information about the batch file. The lion's share of this information
   is that of the batch item. The class is implemented as a singleton so that we don't have to
   deal with yet another global variable.
         
   There are two ways of using this class and it's components. The first involves using try/catch
   blocks around access immediate accesses (as shown below). The second involves wrapping the entire
   program in a try/catch block so that batch file exceptions are processed only in one place.
         
   // Begin by processing the batch file and catching any exceptions generated from that activity.
   // All execptions at this point are fatal, so you will probably want to write a message to the
   // log file and then exit with appropriate exception type.
   try {
       BatchFile::BF()->ProcessBatchFile("MEGA2.BATCH");
   } catch (BatchFileException &e) {
       // Here we give the user a textual description of the exception and exit with the appropriate   value...
       std::cout << "Exception: " << e.what() << std::endl;
       exit(e.getType());
   }
         
   // After processing the batch file you will usually get the derived class associated with a
   // keyword, and possible section (the empty section assumed below) so that you can get information
   // about it. The programmer will need to know that the keyword represents an int. If they are
   // wrong, then a 'std::bad_cast' exception will be thrown. Here the BatchItem information associated
   // with the keyword "Chromosome_Single" is retrieved...
         
   try {
       BatchItemInt& chromosomeSingle = (dynamic_cast<BatchItemInt&>(BatchFile::BF()->getBatchItem("Chromosome_Single")));
   } catch (std::bad_cast) {
       // Only get the exception when casting to a reference.
       // When casting to a pointer you get NULL and no exception.
       std::cout << "Exception casting a B to a C, ooops!" << std::endl;
   }
         
   // It is now possible to ask questions about the batch item...
   int v = chromosomeSingle.getValue();
         
*/
class BatchFile {
 private:
  // Currently the version number is not used...
  double version;
  string batchFileName;
  bool processed;
            
  // VERY IMPORTANT NOTE: A pointer MUST be stored here. If the derived class is stored (which
  // is what will happen) Object "slicing" will occur and the derived class specific members
  // of the object being stored will get sliced off, thus the object stored in the vector just
  // acts a Base class object.
  //
  // A vector of all BatchItem...
  vector<BatchItem *> batchItems;
            
  // The list of sections is used for output mostly so that we can iterate over items by
  // section...
  vector<string> sections;
            
  void initBatchItems();
            
  static BatchFile *bfInstance; // Holds the only instance (singleton) that will exist...
            
  // legacy variables...
  bool batchInputFiles;
  bool batchReorder;
  bool batchTraitSelection;
  bool batchErrorParameters;
            
  // private because this is a singleton...
  BatchFile() {
    version = -1;
    batchFileName = "";
    processed = false;
                
    // Populate the batchItems vector....
    batchItems.clear();
    sections.clear();
    initBatchItems();

    batchInputFiles = false;
    batchReorder = false;
    batchTraitSelection = false;
    batchErrorParameters = false;
  }

 protected:
  // This is how you add a BatchItem to a BatchFile.
  // As a side effect it keeps track of sections so that they can be searched through later...
  void addBatchItem(BatchItem *bi) {
    batchItems.push_back(bi);
    if (find(sections.begin(), sections.end(), bi->getSection()) == sections.end()) sections.push_back(bi->getSection());
  }
            
  // Parse [read and process] the batch file...
  void parse(const string fileName);
            
  // Perform consistency checks on the information that has been read from the batch file...
  void consistencyChecks();
            
  void dumpData () {
    cout << "Size batchItems after creation: " <<  batchItems.size() << "\n";
    for (vector<BatchItem *>::iterator bi_it = batchItems.begin(); bi_it != batchItems.end(); ++bi_it) {
      if ((*bi_it)->getKeyword() == "Error_Except_Loci") cout << "@ dumpData Found Error_Except_Loci: " << *bi_it << "\n";
    }
  }


 public:
  // this implements a singleton since we can only have one batch file per run...
  static BatchFile *BF() {
    if (bfInstance == NULL) bfInstance = new BatchFile();
    return bfInstance;
  }
            
  // used to associate the object with a file.
  // Throws an exception if a file has already been opened or it cannot be opened...
  void ProcessBatchFile(const string fileName) {
    parse(fileName);
    consistencyChecks(); // This actually has some side effects as well
  }
  void CopyToMega2BatchItems();
            
  void setBatchInputFiles() { batchInputFiles = true; }
  void setBatchReorder() { batchReorder = true; }
  void setBatchTraitSelection() { batchTraitSelection = true; }
  void setBatchErrorParameters() { batchErrorParameters = true; }
            
  void setBatchInputFiles(bool v) { batchInputFiles = v; }
  void setBatchReorder(bool v) { batchReorder = v; }
  void setBatchTraitSelection(bool v) { batchTraitSelection = v; }
  void setBatchErrorParameters(bool v) { batchErrorParameters = v; }
            
  ~BatchFile() {
    for (vector<BatchItem *>::iterator it = batchItems.begin(); it != batchItems.end(); it++)
      delete(*it);
  }
            
  bool batchFileOpen() { return processed; }
            
  bool getBatchInputFiles() { return batchInputFiles; }
  bool getBatchReorder() { return batchReorder; }
  bool getBatchTraitSelection() { return batchTraitSelection; }
  bool getBatchErrorParameters() { return batchErrorParameters; }
            
  //
  // Return the BatchItem that is associated with the keyword or keyword/section string...
  // In general, this is how Mega2 will get a BatchItem.
  BatchItem& getBatchItem(const string keyword);
  BatchItem& getBatchItem(const string keyword, const string section);
            
  // Determine if the item has been read from the batch file...
  bool itemRead(const string keyword);
  bool itemRead(const string keyword, const string section);
            
  // Overwrite any existing batch file...
  void overWrite(const string fileName);
            
  // Append to an existing batch file those values that were replaced..
  void appendReplaced(const string fileName);
};

#endif /* BATCH_H */
