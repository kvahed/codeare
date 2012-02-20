template <class T> bool 
resetest (Connector<T>* rc) {

    // OUT:                                                                                                                                                     
    Matrix<cxfl>   meas; // measurement                                                                                                                         
    Matrix<cxfl>   mask; // measurement                                                                                                                         

    // IN:                                                                                                                                                      
    Matrix<cxfl>   txm;  // Transmit maps                                                                                                                       
    Matrix<cxfl>   rxm;  // Receive maps                                                                                                                        

    Matrix<double> b0;   // B0 map                                                                                                                              
    Matrix<double> snro; // SNR optimal image                                                                                                                   
    Matrix<double> bet;  // GRE image mask                                                                                                                      

    Matrix<short> bets;  // BET mask                                                                                                                            
    // Read configuration file and initialise backend ------------                                                                                              

    std::string    cf  = std::string (base + std::string (config));
    rc->ReadConfig (cf.c_str());

    int use_bet;
    rc->Attribute ("use_bet", &use_bet);

    stringstream ss;
    string mef, maf;

    ss << base << rc->Attribute("meas");
    mef = ss.str();
    ss.str("");
    ss << base << rc->Attribute("mask");
    maf = ss.str();

    rc->Init(test);
    // -----------------------------------------------------------                                                                                              

    // Read binary data and transmit to backend ------------------                                                                                              

	IO::RAWRead (meas, mef, std::string("VB15"));
    if (use_bet==1)
		IO::RAWRead (mask, maf, std::string("VB15"));

    rc->SetMatrix ("meas", meas);
    rc->SetMatrix ("mask", mask);
    // ----------------------------------------------------------- 

	// Process data on backend -----------------------------------                                                                                              

    rc->Process(test);
    // -----------------------------------------------------------                                                                                              

    // Get back reconstructed data from backend ------------------                                                                                              

    rc->GetMatrix ("txm",  txm);
    rc->GetMatrix ("rxm",  rxm);
    rc->GetMatrix ("mask", mask);
    rc->GetMatrix ("snro", snro);
    rc->GetMatrix ("b0",   b0);
    rc->GetMatrix ("bets", bets);
    // -----------------------------------------------------------                                                                                              

    // Clear RAM and hangup --------------------------------------                                                                                              

    rc->Finalise(test);
    // -----------------------------------------------------------                                                                                              

    // Write data to a single matlab disk ------------------------                                                                                              

    std::string fname = std::string (base + std::string ("maps.mat"));

#ifdef HAVE_MAT_H
    MATFile* mf = matOpen (fname.c_str(), "w");

    if (mf == NULL) {
        printf ("Error creating file %s\n", fname.c_str());
        return false;
    }

	IO::MXDump  (txm, mf,  "txm", "");
    IO::MXDump  (rxm, mf,  "rxm", "");
    IO::MXDump (snro, mf, "snro", "");
    IO::MXDump   (b0, mf,   "b0", "");
    IO::MXDump (bets, mf, "bets", "");
    IO::MXDump (mask, mf, "mask", "");

    if (matClose(mf) != 0) {
        printf ("Error closing file %s\n",fname.c_str());
        return false;
    }
#endif
    // -----------------------------------------------------------                                                                                              

    return true;

}
