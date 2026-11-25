#include <limits>


template<typename DATATYPE>
DATATYPE load_field(string FIELD, std::ifstream & file)
{
    expected_section(FIELD, file);

    DATATYPE value;

    if(!(file>>value))
    throw simple_error("Invalid field "+FIELD);

    return value;
}



template<typename DATATYPE>
void build_vector(
    string FIELD,
    vector<DATATYPE> & v,
    std::ifstream & file
)
{
    expected_section(FIELD, file);

    int amount = load_int("amount", file);

    if(amount<0)
    throw simple_error("Invalid vector size");

    v.resize(amount);

    for(int i=0; i<amount; i++)
    if(!(file>>v[i]))
    throw simple_error("Invalid vector "+FIELD);
}



template<typename DATATYPE>
void build_vector_from_strings(
    string FIELD,
    vector<DATATYPE> & v,
    std::ifstream & file
)
{
    expected_section(FIELD, file);

    int amount = load_int("amount", file);

    if(amount<0)
    throw simple_error("Invalid vector size");

    v.clear();
    v.reserve(amount);

    file.ignore(
        std::numeric_limits<std::streamsize>::max(),
        '\n'
    );

    string s;

    for(int i=0; i<amount; i++)
    {
        if(!std::getline(file, s))
        throw simple_error("Invalid vector "+FIELD);

        v.emplace_back(s);
    }
}



template<typename DATATYPE>
vector<DATATYPE> load_vector(
    string FIELD,
    int amount,
    std::ifstream & file
)
{
    expected_section(FIELD, file);

    if(amount<0)
    throw simple_error("Invalid vector size");

    vector<DATATYPE> v(amount);

    for(int i=0; i<amount; i++)
    if(!(file>>v[i]))
    throw simple_error("Invalid vector "+FIELD);

    return v;
}



template<typename DATATYPE>
pair<vector<DATATYPE>, vector<DATATYPE>> load_pair_of_vectors(
    string FIELD,
    int amount,
    std::ifstream & file
)
{
    expected_section(FIELD, file);

    if(amount<0)
    throw simple_error("Invalid vector size");

    vector<DATATYPE> v1(amount);
    vector<DATATYPE> v2(amount);

    for(int i=0; i<amount; i++)
    if(!(file>>v1[i]))
    throw simple_error("Invalid first vector "+FIELD);

    for(int i=0; i<amount; i++)
    if(!(file>>v2[i]))
    throw simple_error("Invalid second vector "+FIELD);

    return {std::move(v1), std::move(v2)};
}