#ifndef _ARGUMENTS_MANAGER_H_
#define _ARGUMENTS_MANAGER_H_

using std::map;
using std::string;

class ArgumentsManager
{
protected:
    map<string, string> argumentsValues;
    virtual bool strStartsWith(string s, const char* beggining);
public:
    ArgumentsManager(int argc, char** argv);
    const char* getParam(const char* shortName, const char* longName);
};

#endif
