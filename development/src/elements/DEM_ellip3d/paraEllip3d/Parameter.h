/*
In general, built-in C++ types (ints, floats, characters, etc.) can be transmitted over MPI directly.

For types defined by the standard library (such as std::string or std::vector) and some types in Boost 
(such as boost::variant), the Boost.Serialization library already contains all of the required 
serialization code. In these cases, you need only include the appropriate header from the boost/serialization 
directory.

For types that do not already have a serialization header, you will first need to implement serialization 
code before the types can be transmitted using Boost.MPI.
*/

#ifndef PARAMETER_H
#define PARAMETER_H

#include "realtypes.h"
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <boost/mpi.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

namespace dem { 

  class Parameter {

  public:

    // static function is part of the class, not part of the object, so it is used 
    // like Parameter::getSingleton(). it can only access static members.
    // it is public so that it can be called by others like dem::Parameter::getSingleton() 
    // and it can implicitly call private constructor Parameter().
    static Parameter& getSingleton() {
      static Parameter instance; // instantiated on first use, guaranteed to be destroyed
      return instance;
    }
    
    void readIn(const char *input);
    void writeOut();

  private:
    // constructor must be private to avoid instantiation by others because singleton 
    // should only be instantiated by itself
    Parameter() {} 
    ~Parameter() {}
    // make sure these two are unaccessable to avoid copies of singelton
    Parameter(Parameter const&); // don't implement
    void operator=(Parameter const&); // don't implement

  public:
    std::map<std::string, REAL> parameter;
    std::vector<std::pair<REAL, REAL> > gradation;
    std::map<std::string, std::string> datafile; 
    std::vector<REAL> sigmaPath;
    std::vector<std::size_t> cfdPrintPtcls;

  private:
    friend class boost::serialization::access;
    template<class Archive>
      void serialize(Archive & ar, const unsigned int version) {
      ar & parameter;
      ar & gradation;
      ar & sigmaPath;
      ar & cfdPrintPtcls;
    } 
   
  };

}
#endif
