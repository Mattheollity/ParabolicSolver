#include "basic_io.h"

#include <iomanip>
#include <limits>


template<typename ARRAY_TRIPLETS>
void import_spmat(string FILENAME, ARRAY_TRIPLETS & trp)
{
    auto file = load_file(FILENAME);

    long long amount;

    if(!(file>>amount) || amount<0)
    throw simple_error("Invalid triplet amount");

    trp.clear();

    for(long long i=0; i<amount; i++)
    {
        int row, col;
        double value;

        if(!(file>>row>>col>>value))
        throw simple_error("Invalid triplet");

        trp.emplace_back(row, col, value);
    }
}



template<typename ARRAY_TRIPLETS>
void export_spmat(string FILENAME, const ARRAY_TRIPLETS & trp)
{
    std::ofstream file(FILENAME);

    if(!file)
    throw simple_error(FILENAME+" cannot be opened");

    file<<std::setprecision(
        std::numeric_limits<double>::max_digits10
    );

    file<<trp.size()<<"\n";

    for(const triplet & t : trp)
    file<<t.row()<<" "<<t.col()<<" "<<t.value()<<"\n";
}