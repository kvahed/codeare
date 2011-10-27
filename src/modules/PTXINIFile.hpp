/*
 *  jrrs Copyright (C) 2007-2010 Kaveh Vahedipour
 *                               Forschungszentrum Juelich, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but 
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 *  02110-1301  USA
 */

RRSModule::error_code
PTXWriteSiemensINIFile (Matrix<cplx>& pt, int dimrf, int dimgr, int nc, int sampint, float max_rf, std::string* fname, std::string* orientation) {
	
	FILE* fp = fopen (fname->c_str(), "wb");

	int nt = pt.Dim(COL);

	size_t ci[8] = {1, 2, 3, 4, 5, 6, 7, 0};

	if (fp == NULL)
		return RRSModule::FILE_ACCESS_FAILED;

	// Preamble ------------------------------------------
	
	fprintf (fp, "[pTXPulse]\n"                           );
	fprintf (fp, "\n"                                     );
	fprintf (fp, "NUsedChannels    = %i\n",             nc);
	fprintf (fp, "DimRF            = %i\n",          dimrf);
	fprintf (fp, "DimGradient      = %i\n",          dimgr);
	fprintf (fp, "MaxAbsRF         = %3.3f\n",       100.0); // scaling for RF amplitude
	fprintf (fp, "InitialPhase     = %i\n",            0  );
	fprintf (fp, "Asymmetry        = %1.1f\n",         0.5);
	fprintf (fp, "\n"                                    );
	fprintf (fp, "PulseName        = %s\n",          "ktp");
	fprintf (fp, "Family           = %s\n",          "pTX");
	fprintf (fp, "Comment          = %s\n",          "pTX");
	fprintf (fp, "NominalFlipAngle = %i\n",           90  );
	fprintf (fp, "Samples          = %i\n",             nt);
	fprintf (fp, "AmplInt          = %i\n",          100  );
	fprintf (fp, "AbsInt           = %i\n",          100  );
	fprintf (fp, "PowerInt         = %i\n",          100  );
	fprintf (fp, "MinSlice         = %i\n",            1  );
	fprintf (fp, "MaxSlice         = %i\n",            1  );
	fprintf (fp, "RefGrad          = %i\n",            1  );
	fprintf (fp, "\n"                                     );
	fprintf (fp, "\n"                                     );
	fprintf (fp, "[Txa_SAR_SECTION]\n"                    );
	fprintf (fp, "\n"                                     );
	fprintf (fp, "IsPulseSARchecked = %i\n",           0  );
	fprintf (fp, "TX1_PALI_LIMIT_PEAK=\"%1.1f\"\n",    0.0);
	fprintf (fp, "TX1_PALI_LIMIT_10SEC=\"%1.1f\"\n",   0.0);
	fprintf (fp, "TX1_PALI_LIMIT_6MIN=\"%1.1f\"\n",    0.0);
	fprintf (fp, "\n"                                     );
	fprintf (fp, "TX2_PALI_LIMIT_PEAK=\"%1.1f\"\n",    0.0);
	fprintf (fp, "TX2_PALI_LIMIT_10SEC=\"%1.1f\"\n",   0.0);
	fprintf (fp, "TX2_PALI_LIMIT_6MIN=\"%1.1f\"\n",    0.0);
	fprintf (fp, "\n"                                     );
	fprintf (fp, "TX3_PALI_LIMIT_PEAK=\"%1.1f\"\n",    0.0);
	fprintf (fp, "TX3_PALI_LIMIT_10SEC=\"%1.1f\"\n",   0.0);
	fprintf (fp, "TX3_PALI_LIMIT_6MIN=\"%1.1f\"\n",    0.0);
	fprintf (fp, "\n"                                     );
	fprintf (fp, "TX4_PALI_LIMIT_PEAK=\"%1.1f\"\n",    0.0);
	fprintf (fp, "TX4_PALI_LIMIT_10SEC=\"%1.1f\"\n",   0.0);
	fprintf (fp, "TX4_PALI_LIMIT_6MIN=\"%1.1f\"\n",    0.0);
	fprintf (fp, "\n"                                     );
	fprintf (fp, "TX5_PALI_LIMIT_PEAK=\"%1.1f\"\n",    0.0);
	fprintf (fp, "TX5_PALI_LIMIT_10SEC=\"%1.1f\"\n",   0.0);
	fprintf (fp, "TX5_PALI_LIMIT_6MIN=\"%1.1f\"\n",    0.0);
	fprintf (fp, "\n"                                     );
	fprintf (fp, "TX6_PALI_LIMIT_PEAK=\"%1.1f\"\n",    0.0);
	fprintf (fp, "TX6_PALI_LIMIT_10SEC=\"%1.1f\"\n",   0.0);
	fprintf (fp, "TX6_PALI_LIMIT_6MIN=\"%1.1f\"\n",    0.0);
	fprintf (fp, "\n"                                     );
	fprintf (fp, "TX7_PALI_LIMIT_PEAK=\"%1.1f\"\n",    0.0);
	fprintf (fp, "TX7_PALI_LIMIT_10SEC=\"%1.1f\"\n",   0.0);
	fprintf (fp, "TX7_PALI_LIMIT_6MIN=\"%1.1f\"\n",    0.0);
	fprintf (fp, "\n"                                     );
	fprintf (fp, "TX8_PALI_LIMIT_PEAK=\"%1.1f\"\n",    0.0);
	fprintf (fp, "TX8_PALI_LIMIT_10SEC=\"%1.1f\"\n",   0.0);
	fprintf (fp, "TX8_PALI_LIMIT_6MIN=\"%1.1f\"\n",    0.0);
	fprintf (fp, "\n"                                     );
	fprintf (fp, "\n"                                     );

	// Gradient section ----------------------------------

	float maxg = 0.0;

	for (int i = 0; i < nt; i++)
		for (int j = nc; j < nc+3; j++)
			maxg = (maxg > abs(pt(i,j))) ? maxg : abs(pt(i,j));

	fprintf (fp, "[Gradient]\n"                           );
	fprintf (fp, "\n"                                     );
	
	fprintf (fp, "GradientSamples   = %i\n",            nt);
	fprintf (fp, "MaxAbsGradient[0] = %1.5f %1.5f %1.5f\n", maxg, maxg, maxg);
	fprintf (fp, "\n"                                   );
	
	for (int i = 0; i < nt; i++)
		if (orientation->compare("transversal") == 0 || orientation->compare("t") == 0)
			fprintf (fp, "G[%i]= %.4f	 %.4f	 %.4f \n", i, real(pt(i,nc+0))/maxg, real(pt(i,nc+1))/maxg,  real(pt(i,nc+2))/maxg);
		else if (orientation->compare("sagittal") == 0 || orientation->compare("s") == 0)
			fprintf (fp, "G[%i]= %.4f	 %.4f	 %.4f \n", i, real(pt(i,nc+2))/maxg, real(pt(i,nc+1))/maxg, -real(pt(i,nc+0))/maxg);
	
	for (int j = 0; j < nc; j++) {
		
		fprintf (fp, "\n");
		fprintf (fp, "\n");
		fprintf (fp, "[pTXPulse_ch%i]\n", (int)ci[j]);
		fprintf (fp, "\n");
		
		for (int i = 0; i < nt; i++)
			fprintf (fp, "RF[%i]= %.5f	 %.5f\n", i, abs(pt(i,j))*100.0, (arg(pt(i,j)) >= 0.0) ? arg(pt(i,j)) : 6.28318 + arg(pt(i,j)));
		
	}
    
	fclose (fp);
	
	return RRSModule::OK;
	
}
