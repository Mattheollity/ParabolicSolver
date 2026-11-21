#include "basic_io.h"

std::ifstream load_file(string FILE)
{
    if(!fs::exists(FILE))
    throw simple_error(FILE+" file is not exist");

    return std::ifstream(FILE);
}

void expected_section(string FIELD, std::ifstream & file)
{
    string buf;
    file>>buf;
    if(buf!=FIELD) simple_error(string("Expected section ")+FIELD);
}

template<typename DATATYPE>
DATATYPE load_field(string FIELD, std::ifstream & file)
{
    DATATYPE dtbuf;
    string buf;
    file>>buf;
    if(buf!=FIELD) simple_error(string("Expected field ")+FIELD);
    file>>dtbuf;
    return dtbuf; 
}

int load_int(string FIELD, std::ifstream & file)
{
    int ibuf;
    string buf;
    file>>buf;
    if(buf!=FIELD) simple_error(string("Expected int field ")+FIELD);
    file>>ibuf;
    return ibuf;
}

bool load_bool(string FIELD, std::ifstream & file)
{
    bool bbuf;
    string buf;
    file>>buf;
    if(buf!=FIELD) simple_error(string("Expected bool field ")+FIELD);
    file>>bbuf;
    return bbuf;
}

string load_string(string FIELD, std::ifstream & file)
{
    string buf;
    file>>buf;
    if(buf!=FIELD) simple_error(string("Expected string field ")+FIELD);
    file>>buf;
    return buf;
}

template<typename DATATYPE>
vector<DATATYPE> load_vector(string FIELD, int amount, std::ifstream & file)
{
    vector<DATATYPE> vbuf;
    string buf;
    int ibuf;
    file>>buf;
    if(buf!=FIELD) simple_error(string("Expected vector field ")+FIELD+"with "+to_string(amount)+" elements");
    vbuf.resize(amount);
    for(int i=0;i<amount;i++)
    file>>vbuf[i];
    
    return vbuf;
}

template<typename DATATYPE>
pair<vector<DATATYPE>, vector<DATATYPE>> load_pair_of_vectors(string FIELD, int amount, std::ifstream & file)
{
    vector<DATATYPE> vbuf1;
    vector<DATATYPE> vbuf2;
    string buf;
    int ibuf;
    file>>buf;
    if(buf!=FIELD) simple_error(string("Expected pair of vectors field ")+FIELD+"with "+to_string(amount)+" elements each");
    vbuf1.resize(amount);
    vbuf2.resize(amount);
    for(int i=0;i<amount;i++)
    file>>vbuf1[i];
    for(int i=0;i<amount;i++)
    file>>vbuf2[i];
    
    return {vbuf1, vbuf2};
}

template<typename DATATYPE>
void build_vector(string FIELD, vector<DATATYPE> & v, std::ifstream & file)
{
    expected_section(FIELD, file);

    int amount = load_int("amount", file);
    v.resize(amount);

    for(int i=0;i<amount;i++)
    file>>v[i];
}

template<typename DATATYPE>
void build_vector_from_strings(string FIELD, vector<DATATYPE> & v, std::ifstream & file)
{
    expected_section(FIELD, file);

    string s;
    int amount = load_int("amount", file);
    v.resize(amount);

    for(int i=0;i<amount;i++)
    {
        getline(file, s);
        v[i] = DATATYPE(s);
    }
}

bool file_equality(string FILE1, string FILE2) {
    if (!fs::exists(FILE1)) {
        std::cerr << "Не найден файл: " << FILE1 << '\n';
        return false;
    }
    if (!fs::exists(FILE2)) {
        std::cerr << "Не найден файл: " << FILE2 << '\n';
        return false;
    }

    std::ifstream fa(FILE1, std::ios::binary);
    std::ifstream fb(FILE2, std::ios::binary);
    if (!fa || !fb) return false;

    return std::equal(std::istreambuf_iterator<char>(fa),
                      std::istreambuf_iterator<char>(),
                      std::istreambuf_iterator<char>(fb));
}

void write_matrix(const spmat & m, string s)
{
    std::ofstream file(s+".txt");
    file << std::setprecision(std::numeric_limits<double>::max_digits10);
    for(int i=0;i<m.rows();i++)
    {
        for(int j=0;j<m.cols();j++)
        if(m.coeff(i, j)!=0) file<<m.coeff(i, j)<<" ";
        else file<<" . ";
        
        file<<"\n";
    }
}

void write_matrix_struct(const spmat & m, string s)
{
    std::ofstream file(s+".txt");
    file << std::setprecision(std::numeric_limits<double>::max_digits10);
    for(int i=0;i<m.rows();i++)
    {
        for(int j=0;j<m.cols();j++)
        if(m.coeff(i, j)!=0) file<<"!"<<" ";
        else file<<" . ";
        
        file<<"\n";
    }
}

void write_matrix(mat m, string s)
{
    std::ofstream file(s+".txt");
    file << std::setprecision(std::numeric_limits<double>::max_digits10);

    for(int i=0;i<m.rows();i++)
    {for(int j=0;j<m.cols();j++) file<<m.coeff(i, j)<<" "; file<<"\n";}
}

void write_matrix_struct(mat m, string s)
{
    std::ofstream file(s+".txt");
    file << std::setprecision(std::numeric_limits<double>::max_digits10);

    for(int i=0;i<m.rows();i++)
    {
        for(int j=0;j<m.cols();j++)
        {
            if(m.coeffRef(i,j)!=0) file<<"!"<<" ";
            else file<<"."<<" ";
        }
        file<<"\n";
    }
}
void write_matrix(vect m, string s)
{
    std::ofstream file(s+".txt");
    file << std::setprecision(std::numeric_limits<double>::max_digits10);

    for(int i=0;i<m.rows();i++)
    {for(int j=0;j<m.cols();j++) file<<m.coeff(i, j)<<" "; file<<"\n";}
}

void import_spmat(string FILENAME, spmat & mat)
{
    std::ifstream in(FILENAME);

    long long rows, cols, nnz;
    in >> rows >> cols >> nnz;

    std::list<Eigen::Triplet<double>> triplets;

    for (long long i = 0; i < nnz; ++i) {
        long long row, col;
        double val;
        in >> row >> col >> val;

        triplets.push_back(
            triplet(
                static_cast<int>(row),
                static_cast<int>(col),
                val
            ));
    }

    mat = spmat(static_cast<int>(rows), static_cast<int>(cols));
    mat.setFromTriplets(triplets.begin(), triplets.end());
}

void export_spmat(string FILENAME, const spmat & mat)
{
    std::ofstream out(FILENAME);
    out.precision(17);
    out << mat.rows() << " " << mat.cols() << " " << mat.nonZeros() << "\n";

    for (int j = 0; j < mat.outerSize(); ++j)
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat, j); it; ++it)
    out << it.row() << " " << it.col() << " " << it.value() << "\n";
}