#include <string>
#include <map>
#include <gslpp_complex.h>


/**
 * @class Parameters
 * @brief Parameters Class, used to store the parameters for a given model
 */
class Parameters {
  public:
    static const int NumberOfTypes = 4;
    static const std::string TypeList[NumberOfTypes]; 

    /**
     * @brief Void constructor
     */
    Parameters() {};
    /**
     * Parameters copy constructor
     * @param P parameters to be copied
     */
    Parameters(Parameters& P);
    virtual ~Parameters() {};

    /**
    * get the Integers map
    * @return the Integers map
    */
    std::map<std::string, int> getInts();
    /**
     * get the Doubles map
    * @return the Doubles map
    */
    std::map<std::string, double> getDoubles();
    /**
     * get the Complexes map
    * @return the Complexes map
    */
    std::map<std::string, gslpp::complex> getComplexes();
    /**
     * get the Strings map
    * @return the Strings map
    */
    std::map<std::string, std::string> getStrings();

    /**
     * Set the value of an integer
     * @param s the key
     * @param i the value
     */
    void Set(std::string s, int i);
     /**
     * Set the value of a double
     * @param s the key
     * @param d the value
     */
    void Set(std::string s, double d);
    /**
     * Set the value of a complex
     * @param s the key
     * @param z the value
     */
    void Set(std::string s, gslpp::complex z);
    /**
     * Set the value of a string
     * @param s the key
     * @param t the value
     */
    void Set(std::string s, std::string t);

    /**
     * Get the value of an integer
     * @param s the key
     * @param i the value
     */
    void Get(std::string s, int& i);
     /**
     * Get the value of a double
     * @param s the key
     * @param d the value
     */
    void Get(std::string s, double& d);
    /**
     * Get the value of a complex
     * @param s the key
     * @param z the value
     */
    void Get(std::string s, gslpp::complex& z);
    /**
     * Get the value of a string
     * @param s the key
     * @param t the value
     */
    void Get(std::string s, std::string& t);

  private:
    enum MapType {INT, DOUBLE, COMPLEX, STRING};
    std::map<std::string, int> Ints;
    std::map<std::string, double> Doubles;
    std::map<std::string, gslpp::complex> Complexes;
    std::map<std::string, std::string> Strings;

    void InOtherMaps(std::string s, MapType m);
    void InMap(std::string s, MapType m);
};
