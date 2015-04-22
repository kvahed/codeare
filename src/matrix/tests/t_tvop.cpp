#include "Matrix.hpp"
#include "Creators.hpp"
#include "Print.hpp"
#include "IOContext.hpp"
#include "TVOP.hpp"

int test_2d () {
    Matrix<cxfl> A = phantom<cxfl>(256), B, C;
    TVOP<cxfl> tv ;
    B = tv * A;
    C = tv ->* B;
    IOContext f ("tvop_2d.h5", WRITE);
    fwrite (f, A);
    fwrite (f, B);
    fwrite (f, C);
    fclose (f);
	return 0;
}

int test_3d () {
    Matrix<cxfl> A = phantom<cxfl>(128,128,128), B, C;
    TVOP<cxfl> tv ;
    B = tv * A;
    C = tv ->* B;
    IOContext f ("tvop_3d.h5", WRITE);
    fwrite (f, A);
    fwrite (f, B);
    fwrite (f, C);
    fclose (f);
	return 0;
}

int main (int narg, char** argv) {
    test_2d();
    test_3d();
    return 0;
}
