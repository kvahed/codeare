#include "Workspace.hpp"
#include "Algos.hpp"
#include "Creators.hpp"
#include "Print.hpp"

#define VERBOSE

//using namespace codeare::matrix::io;

template <class T>
bool check_workspace () {
    Matrix<T> A = rand<T>(3,4);
    wspace.SetMatrix ("A", A);
    std::cout << A << std::endl;
    std::cout << wspace.Get<T>("A") << std::endl;
    wspace.Free("A");
    return true;
}

int main (int args, char** argv) {

    check_workspace<cxfl>();
    check_workspace<cxdb>();
    //check_workspace<float>();
    //check_workspace<double>();
    //check_workspace<short>();
    //check_workspace<long>();
    return 0;

}
