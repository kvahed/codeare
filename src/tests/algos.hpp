template <class T> bool
algotest (Connector<T>* rc) {

	Matrix<double> A;
	Matrix<cxdb>   B;
	Matrix<cxfl>   C;
	
	MXRead (A, std::string (base + std::string("real.mat")), "A");
	MXRead (B, std::string (base + std::string("complex.mat")), "B");

	Matrix<size_t> sz;

	sz = size(A);

	std::cout << "  size(A): \n    " << sz ;
	std::cout << "  numel(B): " << numel(B) << std::endl;
	
	C = (Matrix<cxfl>) Sum (B,4);
	sz = size(C);
	std::cout << "  size(C): \n    " << sz ;
	std::cout << "  numel(C): " << numel(C) << std::endl;

	MXDump (C, std::string (base + std::string("sumb3.mat")), "sumb3");

	return true;
	
}
