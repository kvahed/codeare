#include "Matrix.hpp"
#include "Creators.hpp"
#include "Print.hpp"
#include "TVOP.hpp"

int test_2d () {
    Matrix<cxfl> A = phantom<cxfl>(256), B, C;
    TVOP<cxfl> tv ;
    B = tv * A;
    C = tv ->* B;
	return 0;
}

int test_3d () {
    Matrix<cxfl> A = phantom<cxfl>(128,128,128), B, C;
    TVOP<cxfl> tv ;
    B = tv * A;
    C = tv ->* B;
	return 0;
}

int test_5d4 () {
    Matrix<cxfl> A = randn<cxfl>(4,5,3,4,3), B, C;
    TVOP<cxfl> tv = TVOP<cxfl>(0,0,0,1,0) ;
    B = tv * A;
    C = tv ->* B;
	return 0;
}

int test_5d5 () {
    Matrix<cxfl> A = randn<cxfl>(4,5,3,4,3), B, C;
    TVOP<cxfl> tv = TVOP<cxfl>(0,0,0,0,1) ;
    B = tv * A;
    C = tv ->* B;
	return 0;
}

int main (int narg, char** argv) {
    test_2d();
    test_3d();
    test_5d4();
    test_5d5();
    return 0;
}
