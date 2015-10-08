#include <Matrix.hpp>
#include <Creators.hpp>
#include <Algos.hpp>
#include <Print.hpp>

template<class T> void check () {
    std::cout << diff(phantom<T>(8,8),2) << std::endl;
    std::cout << diff(phantom<T>(8,8),2,1) << std::endl;
}

int main (int args, char** argv) {

    check<float>();
/*    check<double>();
    check<cxfl>();
    check<cxdb>();*/
    
    return 0;
    
}



