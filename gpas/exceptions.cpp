#include "exceptions.h"

using namespace Exceptions;

FileNotFoundException::FileNotFoundException(const char* filename)
{
    const char* messageTemplate = "Given file %s hasn't been found.\n";
    int mlen = strlen(messageTemplate);
    int nlen = strlen(filename);
    char* message = new char[mlen + nlen];
    sprintf(message, messageTemplate, filename);
    
    this->message = message;
}

FileNotFoundException::~FileNotFoundException()
{
    delete[] message;
}

char* FileNotFoundException::getMessage()
{
    return message;
}



IndexOutOfRangeException::IndexOutOfRangeException(const char* msg)
{
    this->message = strdup(msg);
}

IndexOutOfRangeException::~IndexOutOfRangeException()
{
    delete[] message;
}

char* IndexOutOfRangeException::getMessage()
{
    return message;
}

ConversionException::ConversionException(const char* msg)
{
    this->message = strdup(msg);
}

ConversionException::~ConversionException()
{
    delete[] message;
}

char* ConversionException::getMessage()
{
    return message;
}
