#include "basic_io.h"


std::ifstream load_file(string FILE)
{
    std::ifstream file(FILE);

    if(!file)
    throw simple_error(FILE+" cannot be opened");

    return file;
}



void expected_section(string FIELD, std::ifstream & file)
{
    string buf;

    if(!(file>>buf) || buf!=FIELD)
    throw simple_error("Expected section "+FIELD);
}



int load_int(string FIELD, std::ifstream & file)
{
    return load_field<int>(FIELD, file);
}



bool load_bool(string FIELD, std::ifstream & file)
{
    return load_field<bool>(FIELD, file);
}



string load_string(string FIELD, std::ifstream & file)
{
    return load_field<string>(FIELD, file);
}