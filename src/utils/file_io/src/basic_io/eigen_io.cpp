#include "eigen_io.h"

#include <algorithm>
#include <iterator>
#include <iomanip>
#include <limits>


bool file_equality(string FILE1, string FILE2)
{
    std::ifstream f1(FILE1, std::ios::binary);
    std::ifstream f2(FILE2, std::ios::binary);

    if(!f1 || !f2) return false;

    return std::equal(
        std::istreambuf_iterator<char>(f1),
        std::istreambuf_iterator<char>(),
        std::istreambuf_iterator<char>(f2),
        std::istreambuf_iterator<char>()
    );
}



void write_matrix(const spmat & m, string s)
{
    std::ofstream file(s+".txt");

    file<<std::setprecision(
        std::numeric_limits<double>::max_digits10
    );

    for(int i=0; i<m.rows(); i++)
    {
        for(int j=0; j<m.cols(); j++)
        if(m.coeff(i, j)!=0.0)
        file<<m.coeff(i, j)<<" ";
        else
        file<<" . ";

        file<<"\n";
    }
}



void write_matrix(const mat & m, string s)
{
    std::ofstream file(s+".txt");

    file<<std::setprecision(
        std::numeric_limits<double>::max_digits10
    );

    for(int i=0; i<m.rows(); i++)
    {
        for(int j=0; j<m.cols(); j++)
        file<<m.coeff(i, j)<<" ";

        file<<"\n";
    }
}



void write_matrix(const vect & v, string s)
{
    std::ofstream file(s+".txt");

    file<<std::setprecision(
        std::numeric_limits<double>::max_digits10
    );

    for(int i=0; i<v.size(); i++)
    file<<v[i]<<"\n";
}



void write_matrix_struct(const spmat & m, string s)
{
    std::ofstream file(s+".txt");

    for(int i=0; i<m.rows(); i++)
    {
        for(int j=0; j<m.cols(); j++)
        if(m.coeff(i, j)!=0.0)
        file<<"! ";
        else
        file<<" . ";

        file<<"\n";
    }
}



void write_matrix_struct(const mat & m, string s)
{
    std::ofstream file(s+".txt");

    for(int i=0; i<m.rows(); i++)
    {
        for(int j=0; j<m.cols(); j++)
        if(m.coeff(i, j)!=0.0)
        file<<"! ";
        else
        file<<". ";

        file<<"\n";
    }
}



void import_spmat(string FILENAME, spmat & m)
{
    auto file = load_file(FILENAME);

    int rows, cols;
    long long amount;

    if(!(file>>rows>>cols>>amount))
    throw simple_error("Invalid matrix header");

    if(rows<0 || cols<0 || amount<0)
    throw simple_error("Invalid matrix size");

    vector<triplet> trp;
    trp.reserve(amount);

    for(long long i=0; i<amount; i++)
    {
        int row, col;
        double value;

        if(!(file>>row>>col>>value))
        throw simple_error("Invalid matrix element");

        if(row<0 || row>=rows || col<0 || col>=cols)
        throw simple_error("Invalid matrix index");

        trp.emplace_back(row, col, value);
    }

    m.resize(rows, cols);

    m.setFromTriplets(trp.begin(), trp.end());
}



void export_spmat(string FILENAME, const spmat & m)
{
    std::ofstream file(FILENAME);

    if(!file)
    throw simple_error(FILENAME+" cannot be opened");

    file<<std::setprecision(
        std::numeric_limits<double>::max_digits10
    );

    file<<m.rows()<<" "<<m.cols()<<" "<<m.nonZeros()<<"\n";

    for(int j=0; j<m.outerSize(); j++)
    for(spmat::InnerIterator it(m, j); it; ++it)
    file<<it.row()<<" "<<it.col()<<" "<<it.value()<<"\n";
}