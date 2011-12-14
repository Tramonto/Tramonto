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
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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

class MatchToTolerance
{
public:
  MatchToTolerance(double tolerance, const IntOrDouble& value) 
  : tolerance_(tolerance),
    value_(value)
  {}
  bool operator()(const IntOrDouble& lhs) {
    if (lhs.isInt() && value_.isInt()) {
      return lhs.getInt() == value_.getInt();
    } else if (lhs.isDouble() && value_.isDouble()) {
      return abs(lhs.getDouble() - value_.getDouble()) < tolerance_;
    } else {
      return false;
    }
  }
private:
  double tolerance_;
  IntOrDouble value_;
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

/**
   \brief Trim the leading whitespace from a string.
*/
void trim(std::string& iostring) {
  iostring.erase(0, iostring.find_first_not_of(" \b\t\r\n"));
}

/**
   \brief Read the "Key Output Parameters" section from the given file.
*/
std::vector<IntOrDouble> readKeyOutputParameters(std::istream& in)
{
  typedef std::vector<IntOrDouble> ResultType;
  char const * const kSectionTag = "Key Output Parameters:";
  char const * const kEndOfSectionTag = ".";
  
  // Read until we hit a line with "Key Output Parameters"
  std::string search(kSectionTag);
  std::string line;
  while (!in.eof() && line != search) {
    std::getline(in, line);
  }

  if (line != search) {
    throw std::runtime_error("Could not find the Key Output Parameters section.");
  }
  
  ResultType result;
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
   \brief Read the actual output parameters from a file.
*/
std::vector<IntOrDouble> readActualParameters(std::istream& in)
{
  // Read to the last line of the file.
  std::string line;
  while (getline(in, line).peek() && !in.eof());

  // Extract the actual parameters from the last line of the file.
  std::istringstream lineread(line);
  std::vector<IntOrDouble> actualParameters;
  std::copy(
    std::istream_iterator<IntOrDouble>(lineread),
    std::istream_iterator<IntOrDouble>(),
    std::back_inserter(actualParameters));
  return actualParameters;
}

/**
   \brief Match the parameters within a tolerance.
   
   The expected parameters that are matched by the actual parameters are removed.
*/
void matchParameters(
  const std::vector<IntOrDouble>& actual, 
  std::vector<IntOrDouble>& ioExpected,
  double tolerance)
{
  typedef std::vector<IntOrDouble>::const_iterator Iter;
  for (Iter it = actual.begin(); it != actual.end(); ++it)
  {
    ioExpected.erase(
      std::remove_if(ioExpected.begin(), ioExpected.end(), MatchToTolerance(tolerance, *it)),
      ioExpected.end());
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
  out << "Usage: " << exename << " <readme-file> <output-dat> <tolerance>" << std::endl;
}

int main(int argc, char* argv[])
{
  // Validate arguments.
  if (argc != 4) {
    printUsage(argv[0], std::cerr);
    return -1;
  }

  std::string readmeFile(argv[1]);
  std::string outputDat(argv[2]);
  double tolerance;
  std::istringstream readtolerance(argv[3]);
  readtolerance >> tolerance;

  try {
    std::ifstream readmeIn(readmeFile.c_str());
    std::vector<IntOrDouble> expected = readKeyOutputParameters(readmeIn);
    std::ifstream outputIn(outputDat.c_str());
    std::vector<IntOrDouble> actual = readActualParameters(outputIn);
    printValues("Expected", expected.begin(), expected.end());
    printValues("Actual", actual.begin(), actual.end()); 
    matchParameters(actual, expected, tolerance);

    if (!expected.empty()) {
      std::cout << "[FAILED]" << std::endl;
      printValues("Unmatched", expected.begin(), expected.end());
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
