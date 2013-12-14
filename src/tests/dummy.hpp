template <class T> bool 
dummytest (Connector<T>* rc) {

    if (rc->Init (test) != codeare::OK) {
        printf ("Intialising failed ... bailing out!"); 
        return false;
    }

    Matrix<float> mf (16,4,1,7);
    Matrix<double> md (3,2);
    Matrix<cxfl> mcf (2,2,8);
    Matrix<cxdb> mcd (3,4);
    Matrix<short> msi (7,7);
    Matrix<long> mli (1,8,1);

    rc->SetMatrix("mf", mf);
    rc->SetMatrix("md", md);
    rc->SetMatrix("mcf", mcf);
    rc->SetMatrix("mcd", mcd);
    rc->SetMatrix("msi", msi);
    rc->SetMatrix("mli", mli);

    rc->Process(test);

    rc->Finalise(test);
    
    return true;

}
