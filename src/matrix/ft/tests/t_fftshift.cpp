#include "DFT.hpp"
#include "IOContext.hpp"

int main (int narg, char** argv) {
    using namespace codeare::matrix::io;
    Matrix<cxfl> A = randn<cxfl>(3,6,4,5), B, C;
    B = fftshift(A,3);
    C = fft(A,1);
    IOContext f ("t_fftshift.mat", MATLAB, WRITE);
    fwrite (f, A);
    fwrite (f, B);
    fwrite (f, C);
    fclose(f);
    return 0;
}
