#ifndef _EXCEPTIONS_H_
#define	_EXCEPTIONS_H_ 1

#include <stdio.h>
#include <string.h>


namespace Exceptions
{

    class Exception
    {
    public:
        virtual char* getMessage() = 0;
    };

    class FileNotFoundException : public Exception
    {
    private:
        char* message;
    public:
        FileNotFoundException(const char* filename);
        ~FileNotFoundException();
        virtual char* getMessage();
    };

    class IndexOutOfRangeException : public Exception
    {
    private:
        char* message;
    public:
        IndexOutOfRangeException(const char* msg);
        ~IndexOutOfRangeException();
        virtual char* getMessage();
    };

    class ConversionException : public Exception
    {
    private:
        char* message;
    public:
        ConversionException(const char* msg);
        ~ConversionException();
        virtual char* getMessage();
    };

}

#endif

