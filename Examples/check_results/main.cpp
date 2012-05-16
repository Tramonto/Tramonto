/**
   This program checks the results from a Tramonto example run. In
   particular, it opens the README file in each directory, finds the
   Key Output Parameters section and ensures each parameter appears in
   the last run of the associated example.
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <cmath>

/**
   \brief Determine if a numeric string is floating point or integer valued.
 */
bool isDouble(const std::string& str) {
  return str.find(".") != std::string::npos;
}

/**
   \brief Class that represents an integer or double value.
*/
class IntOrDouble
{
public:

  IntOrDouble(int i = 0) 
  : isInt_(true),
    int_(i) {}
  
  IntOrDouble(double d)
  : isInt_(false),
    double_(d) {}
  
  IntOrDouble(const IntOrDouble& other)
  : isInt_(other.isInt_),
    int_(other.int_),
    double_(other.double_) {}
  
  bool isInt() const {
    return isInt_;
  }
  
  bool isDouble() const {
    return !isInt_;
  }

  void set(int i) {
    isInt_ = true;
    int_ = i;
  }

  void set(double d) {
    isInt_ = false;
    double_ = d;
  }
  
  int getInt() const {
    if (!isInt()) {
      throw std::logic_error("Object is not an integer.");
    }
    return int_;
  }

  double getDouble() const {
    if (!isDouble()) {
      throw std::logic_error("Object is not a double.");
    }
    return double_;
  }

private:
  bool isInt_;
  int int_;
  double double_;
};

std::ostream& operator<<(std::ostream& ostr, const IntOrDouble& v) {
  if (v.isInt()) {
    ostr << "(int)" << v.getInt();
  } else {
    ostr << "(double)" << v.getDouble();
  }
  return ostr;
}

std::istream& operator>>(std::istream& istr, IntOrDouble& v) {
  std::string value;
  istr >> value;
  std::istringstream read(value);
  if (isDouble(value)) {
    double d;
    read >> d;
    v.set(d);
  } else {
    int i;
    read >> i;
    v.set(i);
  }
  return istr;
}

bool operator==(const IntOrDouble& lhs, const IntOrDouble& rhs) {
  if (lhs.isInt() && rhs.isInt()) {
    return lhs.getInt() == rhs.getInt();
  } else if (lhs.isDouble() && rhs.isDouble()) {
    return lhs.getDouble() == rhs.getDouble();
  } else {
    return false;
  }
}

class MatchToTolerance
{
public:
  MatchToTolerance(double tolerance, int iterationBound, const IntOrDouble& value) 
  : tolerance_(tolerance),
    iterationBound_(iterationBound),
    value_(value)
  {}
  bool operator()(const IntOrDouble& lhs) {
    if (lhs.isInt() && value_.isInt()) {
      return std::abs(lhs.getInt() - value_.getInt()) <= iterationBound_;
    } else if (lhs.isDouble() && value_.isDouble()) {
      return (std::abs(lhs.getDouble() - value_.getDouble()) < tolerance_);
    } else {
      return false;
    }
  }
private:
  double tolerance_;
  int iterationBound_;
  IntOrDouble value_;
};

/**
   \brief Trim the leading whitespace from a string.
*/
void trim(std::string& iostring) {
  iostring.erase(0, iostring.find_first_not_of(" \b\t\r\n"));
}

std::vector<IntOrDouble> readParameterSection(std::istream& in) {
  char const * const kEndOfSectionTag = ".";

  std::vector<IntOrDouble> result;
  std::string line;
  std::getline(in, line);
  trim(line);
  while(!in.eof() && line != std::string(kEndOfSectionTag)) {
    std::string value = line.substr(line.find_last_of("=") + 1, line.size());
    std::istringstream read(value);
    if (isDouble(value)) {
      double d;
      read >> d;
      result.push_back(IntOrDouble(d));
    } else {
      int i;
      read >> i;
      result.push_back(IntOrDouble(i));
    }
    std::getline(in, line);
    trim(line);
  }
  return result;
}

/**
   \brief Read the "Key Output Parameters" section from the given file.
*/
std::list<std::vector<IntOrDouble> > readKeyOutputParameters(std::istream& in)
{
  std::list<std::vector<IntOrDouble> > result;
  char const * const kSectionTag = "Key Output Parameters";
  
  // Read until we hit a line with "Key Output Parameters"
  std::string search(kSectionTag);
  std::string line;
  while (!in.eof()) {
    if (line.find(search) != std::string::npos) {
      result.push_back(readParameterSection(in));
    }
    std::getline(in, line);
  }
  return result;
}

/**
   \brief Read the actual output parameters from a file.
*/
std::list<std::vector<IntOrDouble> > readActualParameters(std::istream& in)
{
  std::list<std::vector<IntOrDouble> > result;
  // Read past the first line of the file.
  std::string line;
  std::getline(in, line);
  // Extract each parameter set.
  for (std::getline(in, line); !in.eof(); std::getline(in, line)) {
    trim(line);
    if (!line.empty()) {
      std::istringstream lineread(line);
      std::vector<IntOrDouble> actualParameters;
      std::copy(
		std::istream_iterator<IntOrDouble>(lineread),
		std::istream_iterator<IntOrDouble>(),
		std::back_inserter(actualParameters));
      result.push_back(actualParameters);
    }
  }
  return result;
}

/**
   \brief Match the parameters within a tolerance.
   
   The expected parameters that are matched by the actual parameters are removed.
*/
void matchParameters(
  const std::vector<IntOrDouble>& actual, 
  std::vector<IntOrDouble>& ioExpected,
  double tolerance,
  int iterationBound)
{
  typedef std::vector<IntOrDouble>::const_iterator Iter;
  for (Iter it = actual.begin(); it != actual.end(); ++it)
  {
    std::vector<IntOrDouble>::iterator newEnd 
      = std::remove_if(ioExpected.begin(),
		       ioExpected.end(),
		       MatchToTolerance(tolerance, iterationBound, *it));
    ioExpected.erase(newEnd, ioExpected.end());
  }
}

template<typename IterT>
void printValues(const std::string& label, IterT begin, IterT end)
{
  std::cout << label << ": ";
  std::copy(
	    begin, end,
	    std::ostream_iterator<IntOrDouble>(std::cout, ", "));
  std::cout << std::endl;
}

/**
   \brief Print the program usage instructions.
   
   \param exename The name of the executable as invoked.
   \param out The ostream to print to.
 */
void printUsage(const std::string& exename, std::ostream& out)
{
  out << "Usage: " << exename << " <readme-file> <output-dat> <floating-point-tolerance> <iteration-bound>" << std::endl;
}

int main(int argc, char* argv[])
{
  // Validate arguments.
  if (argc != 5) {
    printUsage(argv[0], std::cerr);
    return -1;
  }

  std::string readmeFile(argv[1]);
  std::string outputDat(argv[2]);
  double tolerance;
  std::istringstream readtolerance(argv[3]);
  readtolerance >> tolerance;
  int iterationBound;
  std::istringstream readIterationBound(argv[4]);
  readIterationBound >> iterationBound;

  std::cout << "Expected parameters file is " << readmeFile << std::endl;
  std::cout << "Data file is " << outputDat << std::endl;
  std::cout << "Floating point comparison tolerance is " << tolerance << std::endl;
  std::cout << "Iteration bound is " << iterationBound << std::endl;
  

  try {
    typedef std::list<std::vector<IntOrDouble> > ParameterSet;
    std::ifstream readmeIn(readmeFile.c_str());
    if (!readmeIn.good()) {
      std::cerr << "Failed to read from " << readmeFile << std::endl;
      return 1;
    }
    ParameterSet expected = readKeyOutputParameters(readmeIn);
    if (expected.empty()) {
      std::cerr << "Failed to locate Key Output Parameters section" << std::endl;
      return 1;
    }
    std::ifstream outputIn(outputDat.c_str());
    if (!outputIn.good()) {
      std::cerr << "Failed to read from " << outputDat << std::endl;
      return 1;
    }
    ParameterSet actual = readActualParameters(outputIn);
    if (actual.empty()) {
      std::cerr << "Failed to locate parameters in output file";
    }
    
    std::cout << "Parameter Sets" << std::endl;
    for (ParameterSet::iterator expectedValues = expected.begin(); expectedValues != expected.end(); ++expectedValues) {
      printValues("Expected", expectedValues->begin(), expectedValues->end());
    }
    for (ParameterSet::iterator actualValues = actual.begin(); actualValues != actual.end(); ++actualValues) {
      printValues("Actual", actualValues->begin(), actualValues->end());
    }
    
    // Try to match one of the actual parameter sets to one of the expected parameter sets.
    bool foundMatch = false;
    for (ParameterSet::iterator expectedValues = expected.begin(); expectedValues != expected.end(); ++expectedValues) {
      for (ParameterSet::iterator actualValues = actual.begin(); actualValues != actual.end(); ++actualValues) {
	matchParameters(*actualValues, *expectedValues, tolerance, iterationBound);
	if (expectedValues->empty()) {
	  foundMatch = true;
	}
      }
    }

    if (!foundMatch) {
      std::cout << "[FAILED]" << std::endl;
      for (ParameterSet::iterator expectedValues = expected.begin(); expectedValues != expected.end(); ++expectedValues) {
	std::cout << "Expected parameter set " << std::distance(expected.begin(), expectedValues) << std::endl;
	printValues("  Unmatched", expectedValues->begin(), expectedValues->end());
      }
      return 1;
    } else {
      std::cout << "[PASSED]" << std::endl;
      return 0;
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
}
