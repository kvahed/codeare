#include "Matrix.hpp"
#include "Creators.hpp"
#include "Print.hpp"


int main (int narg, const char** argv) {

	typedef Range<false> R;
	typedef Range<true> CR;

	Matrix<double> M1(10, 1), M2, M3(M1), M4(3, 4), M6(4, 3), M5(3, 4, 2), M7(3, 4, 2), M8, M9, M10, M11;
    M1 = 0.0;
    M2 = M1;

    for (size_t i = 1; i < M1.Size(); ++i)
        M1[i] = i;
    for (size_t i = 1; i < M4.Size(); ++i)
        M4[i] = i;
    for (size_t i = 1; i < M5.Size(); ++i)
        M5[i] = i;

    M3(R(7,-2,2)) = M1(CR(2,4));
    M6(R(3,-2,0),R(0,1)) = M4(CR(0,1),CR(1,2));
    M7(R(2,-2,0),R(0,1),R(0,0)) =  M5(CR(0,1),CR(1,2),CR(0,0));
    M8 = M6(CR(1,-1,0),CR(1,2));
    M6(R(0,1),R(1,2)) = M8;
    M9 = M6(CR("1:-1:0"),CR("1:2"));
    M10 = M6(CR("1:-1:0,1"),CR("1:2"));
    M11 = M6(CR(),CR()) * M6(CR(),CR());
//    M12 = M5(CR("1:end"),CR("1:end"));
//    M11 = M6(CR(),CR()) + M6(CR(),CR());
  //  M11 = M6(CR(),CR()) / M6(CR(),CR());
    M11 = M6 * M6(CR(),CR());
    std::cout << M3 << std::endl<< std::endl;
    std::cout << M4 << std::endl<< std::endl;
    std::cout << M6 << std::endl<< std::endl;
    std::cout << M5 << std::endl<< std::endl;
    std::cout << M7 << std::endl<< std::endl;
    std::cout << M8 << std::endl<< std::endl;
    std::cout << M9 << std::endl<< std::endl;
    std::cout << M10 << std::endl<< std::endl;
    std::cout << M11 << std::endl<< std::endl;
    return 0;
}

