#include <string>
#include <map>

#include "arguments_manager.h"

ArgumentsManager::ArgumentsManager(int argc, char** argv)
{
    string s;
    string val;


    for (int i = 1; i < argc; i++)
    {
        s = argv[i];

        if (strStartsWith(s, "--"))
            s = s.substr(2, s.length() - 2);
        else if (strStartsWith(s, "-"))
            s = s.substr(1, s.length() - 1);
        else
            continue;

        if (i + 1 >= argc)
        {
            argumentsValues[s] = "true";
            continue;
        }

        if (argv[i + 1][0] == '-')
            continue;

        val = argv[i + 1];

        argumentsValues[s] = val;
    }
}

bool ArgumentsManager::strStartsWith(string s, const char* beggining)
{
    char* value = (char*)s.c_str();
    char* tchar = (char*)beggining;
    while (*value)
    {
        if (!*tchar)
            return true;
        if (*value != *tchar)
            return false;
        tchar++;
        value++;
    }
    return false;
}

const char* ArgumentsManager::getParam(const char* shortName, const char* longName)
{
    //string val = argumentsValues[shortName];
    map<string, string>::const_iterator tester = argumentsValues.find(shortName);
    if (tester != argumentsValues.end())
        return tester->second.c_str();
    tester = argumentsValues.find(longName);
    if (tester != argumentsValues.end())
        return tester->second.c_str();
    return NULL;
}
