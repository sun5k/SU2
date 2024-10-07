/*!
 * \file trans_correlations.hpp
 * \brief Numerics class for the LM model's correlation functions.
 * \version 7.5.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

/*!
 * \class TransLMCorrelations
 * \brief Class for LM model's correlation functions.
 * \ingroup SourceDiscr
 * \author A. Rausa.
 */
class TransLMCorrelations {
 private:

  LM_ParsedOptions options;

 public:

  /*!
   * \brief Set LM options.
   * \param[in] val_options - LM options structure.
   */
  void SetOptions(const LM_ParsedOptions val_options){
    options = val_options;
  }

  /*!
   * \brief Compute Re_theta_c from correlations.
   * \param[in] Tu - Turbulence intensity.
   * \param[in] Re_theta_t - Re_theta_t (TransVar[1]).
   * \param[out] rethetac - Corrected value for Re_theta.
   */
  su2double ReThetaC_Correlations(const su2double Tu, const su2double Re_theta_t) const {

    su2double rethetac = 0.0;

    switch (options.Correlation) {
      case TURB_TRANS_CORRELATION::MALAN: {
        rethetac = min(0.615 * Re_theta_t + 61.5, Re_theta_t);
        break;
      }

      case TURB_TRANS_CORRELATION::SULUKSNA: {
        rethetac = min(0.1 * exp(-0.0022 * Re_theta_t + 12), 300.0);
        break;
      }

      case TURB_TRANS_CORRELATION::KRAUSE: {
        rethetac = 0.91 * Re_theta_t + 5.32;
        break;
      }

      case TURB_TRANS_CORRELATION::KRAUSE_HYPER: {
        const su2double FirstTerm = -0.042 * pow(Tu, 3);
        const su2double SecondTerm = 0.4233 * pow(Tu, 2);
        const su2double ThirdTerm = 0.0118 * pow(Tu, 1);
        rethetac = Re_theta_t / (FirstTerm + SecondTerm + ThirdTerm + 1.0744);
        break;
      }

      case TURB_TRANS_CORRELATION::MEDIDA_BAEDER: {
        const su2double FirstTerm = 4.45 * pow(Tu, 3);
        const su2double SecondTerm = 5.7 * pow(Tu, 2);
        const su2double ThirdTerm = 1.37 * pow(Tu, 1);
        rethetac = (FirstTerm - SecondTerm + ThirdTerm + 0.585) * Re_theta_t;
        break;
      }

      case TURB_TRANS_CORRELATION::MEDIDA: {
        rethetac = 0.62 * Re_theta_t;
        break;
      }

      case TURB_TRANS_CORRELATION::MENTER_LANGTRY: {
        if (Re_theta_t <= 1870) {
          const su2double FirstTerm = (-396.035 * pow(10, -2));
          const su2double SecondTerm = (10120.656 * pow(10, -4)) * Re_theta_t;
          const su2double ThirdTerm = (-868.230 * pow(10, -6)) * pow(Re_theta_t, 2);
          const su2double ForthTerm = (696.506 * pow(10, -9)) * pow(Re_theta_t, 3);
          const su2double FifthTerm = (-174.105 * pow(10, -12)) * pow(Re_theta_t, 4);
          rethetac = FirstTerm + SecondTerm + ThirdTerm + ForthTerm + FifthTerm;
        } else {
          rethetac = Re_theta_t - (593.11 + 0.482 * (Re_theta_t - 1870.0));
        }

        break;
      }
      case TURB_TRANS_CORRELATION::DEFAULT:
        SU2_MPI::Error("Transition correlation is set to DEFAULT but no default value has ben set in the code.",
                       CURRENT_FUNCTION);
        break;
    }

    return rethetac;
  }

  /*!
   * \brief Compute FLength from correlations.
   * \param[in] Tu - Turbulence intensity.
   * \param[in] Re_theta_t - Re_theta_t (TransVar[1]).
   * \param[out] F_length1 - Value for the F_length1 variable.
   */
  su2double FLength_Correlations(const su2double Tu, const su2double Re_theta_t) const {
    su2double F_length1 = 0.0;

    switch (options.Correlation) {
      case TURB_TRANS_CORRELATION::MALAN: {
        F_length1 = min(exp(7.168 - 0.01173 * Re_theta_t) + 0.5, 300.0);
        break;
      }

      case TURB_TRANS_CORRELATION::SULUKSNA: {
        const su2double FirstTerm = -pow(0.025 * Re_theta_t, 2) + 1.47 * Re_theta_t - 120.0;
        F_length1 = min(max(FirstTerm, 125.0), Re_theta_t);
        break;
      }

      case TURB_TRANS_CORRELATION::KRAUSE: {
        F_length1 = 3.39 * Re_theta_t + 55.03;
        break;
      }

      case TURB_TRANS_CORRELATION::KRAUSE_HYPER: {
        if (Tu <= 1.) {
          F_length1 = log(Re_theta_t + 1) / Tu;
        } else {
          const su2double FirstTerm = 0.2337 * pow(Tu, 2);
          const su2double SecondTerm = -1.3493 * pow(Tu, 1);
          F_length1 = log(Re_theta_t + 1) * (FirstTerm + SecondTerm + 2.1449);
        }
        break;
      }

      case TURB_TRANS_CORRELATION::MEDIDA_BAEDER: {
        const su2double FirstTerm = 0.171 * pow(Tu, 2);
        const su2double SecondTerm = 0.0083 * pow(Tu, 1);
        F_length1 = (FirstTerm - SecondTerm + 0.0306);
        break;
      }

      case TURB_TRANS_CORRELATION::MEDIDA: {
        F_length1 = 40;
        break;
      }

      case TURB_TRANS_CORRELATION::MENTER_LANGTRY: {
        if (Re_theta_t < 400) {
          F_length1 = 39.8189 + (-119.270 * pow(10, -4)) * Re_theta_t +
                      (-132.567 * pow(10, -6)) * Re_theta_t * Re_theta_t;
        } else if (Re_theta_t < 596) {
          F_length1 = 263.404 + (-123.939 * pow(10, -2)) * Re_theta_t +
                      (194.548 * pow(10, -5)) * pow(Re_theta_t, 2) +
                      (-101.695 * pow(10, -8)) * pow(Re_theta_t, 3);
        } else if (Re_theta_t < 1200) {
          F_length1 = 0.5 - (3.0 * pow(10, -4)) * (Re_theta_t - 596.0);
        } else {
          F_length1 = 0.3188;
        }
        break;
      }
      case TURB_TRANS_CORRELATION::DEFAULT:
        SU2_MPI::Error("Transition correlation is set to DEFAULT but no default value has ben set in the code.",
                       CURRENT_FUNCTION);
        break;
    }

    return F_length1;
  }
};



/*!
 * \class TransAFMTCorrelations
 * \brief Class for AFMT model's correlation functions.
 * \ingroup SourceDiscr
 * \author S. Kang.
 */
class TransAFMTCorrelations {
 private:

  AFMT_ParsedOptions options;

 public:

  /*!
   * \brief Set AFMT options.
   * \param[in] val_options - AFMT options structure.
   */
  void SetOptions(const AFMT_ParsedOptions val_options){
    options = val_options;
  }

  /*!
   * \brief Compute H12 from correlations.
   * \param[in] HL - Local Shape Factor.
   * \param[in] T_over_T0 - T/T0.
   * \param[in] M_e - Edge Mach number.
   * \param[in] Tw_over_Te - Tw/Te.
   * \param[out] H12 - Integrated Shape Factor.
   */
  su2double H12_Correlations(const su2double HL, const su2double T_over_T0, const su2double M_e, const su2double Tw_over_Te) const {

    su2double H12 = 0.0;
    su2double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0;
    su2double b1 = 0.0, b2 = 0.0, b3 = 0.0, b4 = 0.0;

    switch (options.Correlation) {
      case AFMT_CORRELATION::Liu2023: {
        //H12 = 2.816 * Tw_over_Te + 0.1189 * HL + 0.1810 * M_e * M_e - 0.2772;
        H12 = (9.2706e+1 * pow(T_over_T0,2) + 1.2516e+1 * T_over_T0 * M_e - 2.4463e+2 * T_over_T0 - 1.6690e+1 * M_e + 2.0717e+2) * T_over_T0;
        H12 = H12 + 1.7872 * pow(M_e,2) + 5.4920 * M_e - 5.3168e+1;
        H12 = H12 * HL * T_over_T0;
        break;
      }

      case AFMT_CORRELATION::sok: {
        /*
        H12 = (-1.8473 * HL - 1.5789e-1 * T_over_T0 + 3.0071) * pow(M_e,2) + (2.8563 * pow(T_over_T0,2) - 3.4846 * T_over_T0 + 1.9417e+1 * HL - 1.0705e+1 ) * M_e;
        H12 = H12 - 6.259e+1 * pow(HL,2) + 3.6986e+1*HL;
        H12 = H12 + 2.8268e+1 * pow(T_over_T0,3) - 7.5314e+1 * pow(T_over_T0,2) + 6.9668e+1 * T_over_T0;
        H12 = H12 - 1.7524e+1;
        H12 = H12 * HL * T_over_T0;
        */
        H12 = (1.0696e+1 * pow(T_over_T0,2) + 1.4386e+1 * T_over_T0 * M_e - 2.7767e+2 * T_over_T0 - 1.8548e+1 * M_e + 2.2921e+2) * T_over_T0;
        H12 = H12 + 1.7856 * pow(M_e,2) + 5.9571 * M_e - 5.7666e+1;
        H12 = H12 * HL * T_over_T0;
        break;
      }

      case AFMT_CORRELATION::DEFAULT:
        SU2_MPI::Error("Transition correlation is set to DEFAULT but no default value has ben set in the code.",
                       CURRENT_FUNCTION);
        break;
    }
    H12 = max(2.2, H12);
    return H12;
  }

  /*!
   * \brief Compute Hk from correlations.
   * \param[in] HL - Local Shape Factor.
   * \param[in] H12 - Integreated Shape Factor.
   * \param[in] M_e - T/T0.
   * \param[out] Hk - Kinetic Shape Factor.
   */
  su2double Hk_Correlations(const su2double HL, const su2double H12, const su2double M_e) const {
    su2double Hk = 0.0;
    su2double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0;
    su2double b1 = 0.0, b2 = 0.0, b3 = 0.0, b4 = 0.0;

    switch (options.Correlation) {
      case AFMT_CORRELATION::Liu2023: {
        /*
        a1 = 4.957e-3 * pow(M_e,3) - 6.297e-2 * pow(M_e,2) -5.984e-2 * M_e - 2.334;
        a2 = -7.960e-3 * pow(M_e,3) + 1.095e-1 * pow(M_e,2) + 2.831e-3 * M_e + 5.209;
        a3 = 3.327e-3 * pow(M_e,3) - 4.705e-2 * pow(M_e,2) + 1.867e-2 * M_e + 4.361e-1;
        Hk = a1 * pow(HL,2) + a2 * HL + a3;
        */
        a1 = 2.642158e-5 * pow(M_e,5) - 1.047983e-3 * pow(M_e,4) + 1.688650e-2 * pow(M_e,3) - 1.392065e-1 * pow(M_e,2) + 5.936119e-1 * M_e - 1.065709;
        a2 = 1.140608e-4 * pow(M_e,4) - 9.542641e-3 * pow(M_e,3) + 2.129670e-1 * pow(M_e,2) - 1.968956 * M_e + 7.184406;
        a3 = 4.129277e-3 * pow(M_e,4) - 1.256703e-1 * pow(M_e,3) + 1.403225 * pow(M_e,2) - 6.410845 * M_e + 1.309789e+1;
        Hk = a1 * pow(H12,2) + a2 * H12 + a3;
        Hk = log(Hk);
        break;
      }

      case AFMT_CORRELATION::sok: {
        a1 = 2.642158e-5 * pow(M_e,5) - 1.047983e-3 * pow(M_e,4) + 1.688650e-2 * pow(M_e,3) - 1.392065e-1 * pow(M_e,2) + 5.936119e-1 * M_e - 1.065709;
        a2 = 1.140608e-4 * pow(M_e,4) - 9.542641e-3 * pow(M_e,3) + 2.129670e-1 * pow(M_e,2) - 1.968956 * M_e + 7.184406;
        a3 = 4.129277e-3 * pow(M_e,4) - 1.256703e-1 * pow(M_e,3) + 1.403225 * pow(M_e,2) - 6.410845 * M_e + 1.309789e+1;
        Hk = a1 * pow(H12,2) + a2 * H12 + a3;
        Hk = log(Hk);
        break;
      }

      case AFMT_CORRELATION::DEFAULT:
        SU2_MPI::Error("Transition correlation is set to DEFAULT but no default value has ben set in the code.",
                       CURRENT_FUNCTION);
        break;
    }
    /*
    if(Hk == nan){
      cout << "NaN is tected." << endl;
    }
    */
    Hk = min(max(0.0, Hk), 4.0);    
    return Hk;
  }

  /*!
   * \brief Compute Hk from correlations.
   * \param[in] HL - Local Shape Factor.
   * \param[in] H12 - Integreated Shape Factor.
   * \param[in] M_e - T/T0.
   * \param[out] Hk - Kinetic Shape Factor.
   */
  su2double RevRet_Correlations(const su2double H12, const su2double M_e) const {
    su2double RevRet = 0.0;
    su2double a1 = 0.0, a2 = 0.0, a3 = 0.0;

    switch (options.Correlation) {
      case AFMT_CORRELATION::Liu2023: {
        a1 = 0.0008 * pow(M_e,2) + 0.0932 * M_e + 0.1109;
        a2 = -0.0356 * pow(M_e,2) - 0.1249 * M_e + 0.9068;
        a3 = 0.0924 * pow(M_e,2) - 0.7116 * M_e + 2.3833;

        RevRet = a1 * pow(log(H12),2) + a2 * log(H12) + a3;
        break;
      }

      case AFMT_CORRELATION::sok: {
        a1 = 0.0008 * pow(M_e,2) + 0.091 * M_e + 0.1253;
        a2 = -0.0343 * pow(M_e,2) - 0.1322 * M_e + 0.8943;
        a3 = 0.0887 * pow(M_e,2) - 0.6799 * M_e + 2.3288;

        RevRet = a1 * pow(log(H12),2) + a2 * log(H12) + a3;
        break;
      }

      case AFMT_CORRELATION::DEFAULT:
        SU2_MPI::Error("Transition correlation is set to DEFAULT but no default value has ben set in the code.",
                       CURRENT_FUNCTION);
        break;
    }
    /*
    if(RevRet == nan){
      cout << "NaN is tected." << endl;
    }
    */
    RevRet = max(min(20.0, RevRet),0.5);
    return RevRet;
  }





  /*!
   * \brief Compute dN/dRet from correlations.
   * \param[in] H12 - Integreated Shape Factor.   
   * \param[in] M_e - Edge Mach number.
   * \param[out] dNdRet - N factor gradient for Mack 2nd mode.
   */
  su2double dNdRet_Correlations(const su2double H12, const su2double M_e) const {
    su2double dNdRet = 0.0;
    su2double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0;

    switch (options.Correlation) {
      case AFMT_CORRELATION::Liu2023: {
        a1 = -3.030e-3 * pow(M_e,2) - 3.827e-2 * M_e + 7.520e-1;
        a2 = -1.255e-3 * pow(M_e,2) + 1.581e-1 * M_e - 1.605;
        a3 = 2.958e-5 *  pow(M_e,2) - 1.277e-3 * M_e + 1.164e-2;
        a4 = -9.993e-4 * pow(M_e,3) + 1.769e-2 * pow(M_e,2) - 9.183e-2 * M_e + 1.115e-1;
        dNdRet = a1 * exp(a2 * H12) + a3 * exp(a4 * H12);
        dNdRet = min(dNdRet, 0.02);
        break;
      }

      case AFMT_CORRELATION::sok: {
        /*
        a1 = 0.0;
        a2 = 0.0;
        a3 = 0.0;
        dNdRet = a1 * pow(H12,2) + a2 * H12 + a3;
        */
        /* ver_1
        a1 = - 5.2186e-7 * pow(M_e,3) + 1.1325e-5 * pow(M_e,2) - 8.2659e-5 * M_e + 2.0309e-4;
        a2 = + 1.5148e-5 * pow(M_e,3) - 3.5602e-4 * pow(M_e,2) + 2.8322e-3 * M_e - 7.6569e-3;
        a3 = + 3.1036e-4 * pow(M_e,2) - 5.6516e-3 * M_e + 2.9353e-2;
        */
        // ver_2
        /*
        a1 = 0.007406 * exp(-1.331 * M_e);
        a2 = -0.05271 * exp(-0.9117 * M_e);
        a3 = 0.03892 * exp(-0.2769 * M_e);
        dNdRet = a1 * pow(H12, 2) + a2 * pow(H12, 1) + a3;
        */
        // ver2_3 custom
        /*
        a1 = 1.2393e-04 * exp(-4.4071e-01 * M_e);
        a2 = -4.3843e-03 * exp(-3.4754e-01 * M_e);
        a3 = -3.4257e-04 * pow(M_e,2) +3.4118e-03 * M_e + 3.2536e-03;
        dNdRet = a1 * pow(H12, 2) + a2 * pow(H12, 1) + a3;
        */
       //ver2_4 custom2
       /*
        a1 =(0.005211*pow(M_e, 2)-0.0495*pow(M_e,1)+0.1172)/(pow(M_e,2)-10.77*pow(M_e,1)+29.7); 
        a2 = -0.3235*pow(M_e,2)+2.634*pow(M_e,1)-4.992;
        if(M_e > 5.5 ){
          a2 = 0.00639*pow(M_e,4)-0.08947*pow(M_e,3)-0.02045*pow(M_e,2)+4.613*pow(M_e,1)-16.01;
        }
         
        a3 = -0.002592*pow(M_e,3)+0.05448*pow(M_e,2)-0.3797*pow(M_e,1)+0.8783;
        a4 = -0.04386*pow(M_e,4)+0.9824*pow(M_e,3)-8.222*pow(M_e,2)+30.61*pow(M_e,1)-43.02;
        if(M_e > 6.5){
          a4 = 0.195*pow(M_e,3)-4.201*pow(M_e,2)+30.05*pow(M_e,1)-71.37;
        }
        dNdRet = a1 * exp(a2 * H12) + a3 * exp(a4 * H12);
        */
        //ver 2.5
        /*
        a1 = -0.000128 * pow(M_e,3)+0.002182 *pow(M_e,2)-0.012479*pow(M_e,1)+ 0.024005;
        a2 = 0.0038133 * pow(M_e,3)- 0.064 *pow(M_e,2)+ 0.36213*pow(M_e,1)- 0.69354;
        a3 = -0.0026667 * pow(M_e,3) + 0.052 *pow(M_e,2)- 0.34333*pow(M_e,1)+0.84;
        a4 = -3.6 *pow(M_e, 3) + 55.6 *pow(M_e, 2)- 278.1*pow(M_e,1)+453.0;
        if(M_e > 6.0) {
          a1 = 0.0000026667 * pow(M_e,4) - 0.000073333 * pow(M_e,3) + 0.00074533 * pow(M_e,2) - 0.0033227 * pow(M_e,1) + 0.005523;
          a2 = 0.0008 * pow(M_e, 4) - 0.021733*pow(M_e,3) + 0.2184*pow(M_e,2) - 0.96277*pow(M_e,1) + 1.5707;
          a3 = 0.022667 * pow(M_e,4) - 0.584*pow(M_e,3) + 5.6203*pow(M_e,2) - 23.919*pow(M_e,1) + 38.026;
          a4 = 2.8667 *pow(M_e, 4) - 77.8*pow(M_e, 3) + 790.98*pow(M_e, 2) - 3568.9*pow(M_e, 1) + 6036.0;
        }
        if(pow(H12,2) == - a4){
          dNdRet = 0.0;
        }
        else {
          dNdRet = (a1 * pow(H12,2) + a2 * pow(H12,1) + a3) / ( pow(H12,2) + a4 );
        }
        */

        //ver 2.6
        /*

        a1 =  1.73466667E-01*pow(M_e,3) - 2.90080000E+00*pow(M_e,2) + 1.62186333E+01*pow(M_e,1) - 3.02698000E+01;
        a2 = -1.44133333E+00*pow(M_e,3) + 2.52080000E+01*pow(M_e,2) - 1.46579667E+02*pow(M_e,1) + 2.83568000E+02;
        a3 = -1.18933333E+01*pow(M_e,3) + 1.99560000E+02*pow(M_e,2) - 1.08692667E+03*pow(M_e,1) + 1.92336000E+03;
        a4 = -1.10000000E+01*pow(M_e,3) + 2.65600000E+02*pow(M_e,2) - 2.00705000E+03*pow(M_e,1) + 4.85570000E+03;
        if(M_e > 6.0) {
          a1 = 6.88400000E-02*pow(M_e,4) - 1.83552000E+00*pow(M_e,3) + 1.80489100E+01*pow(M_e,2) - 7.75947800E+01*pow(M_e,1) + 1.23145600E+02;
			    a2 = 4.03333333E+00*pow(M_e,4) - 1.08833333E+02*pow(M_e,3) + 1.10169167E+03*pow(M_e,2) - 4.95464167E+03*pow(M_e,1) + 8.34800000E+03;
			    a3 = -3.54550667E+01*pow(M_e,4) + 1.00140853E+03*pow(M_e,3) - 1.05201374E+04*pow(M_e,2) + 4.87920545E+04*pow(M_e,1) - 8.43648560E+04;
			    a4 = 1.21080000E+01*pow(M_e,4) - 3.51690667E+02*pow(M_e,3) + 3.80661700E+03*pow(M_e,2) - 1.81481193E+04*pow(M_e,1) + 3.21227200E+04;

        }
         
        if(pow(H12,2) + a3 * pow(H12,1) == - a4){
          dNdRet = 0.000000000001;
        }
        else {
          dNdRet = (a1 * pow(H12,1) + a2 ) / ( pow(H12,2) + a3 * pow(H12,1) + a4 );
        }
        */
        //ver 3.0
        /*
        a1 = 0.007406 * exp(-1.331 * M_e);
        a2 = -0.05271 * exp(-0.9117 * M_e);
        a3 = 0.03892 * exp(-0.2769 * M_e);
        dNdRet = a1 * pow(H12, 2) + a2 * pow(H12, 1) + a3;
        */
        //ver 3.1
        /*
        a1 = 1.6251E-02 * exp(-1.1152E+00 * M_e) + 4.1051E-06 * exp(-1.0866E-02 * M_e) ;
        a2 = -2.5929E-04 * pow(M_e,2) + 4.0741E-03*pow(M_e,1) -1.6397E-02;
        a3 = 4.0055E-04 * pow(M_e,3) -7.2034E-03 *pow(M_e,2)+3.7555E-02*pow(M_e,1)-3.3531E-02;
        dNdRet = a1 * pow(H12, 2) + a2 * pow(H12, 1) + a3;
        */
        //ver 3.2
        /*
        a1 = 4.1481E-03* pow(M_e,3) - 6.6463E-02* pow(M_e,2) + 3.2469E-01* pow(M_e,1) - 4.3057E-01;
        a2 = -3.5286E-02* pow(M_e,3) + 5.5082E-01* pow(M_e,2) - 2.7029E+00* pow(M_e,1) + 3.8855E+00;
        a3 = 0.0;
        if(M_e >6.0){
          a3 = -1.9143E-03* pow(M_e,2) + 2.8340E-02* pow(M_e,1) - 1.0118E-01;
        }
        dNdRet = a1 * exp(a2 * H12) + a3 ;
        */
        //ver 3.3
        /*
        a1 = -2.266667E-01 * pow(M_e,3) +3.020000E+00* pow(M_e,2) -1.325333E+01* pow(M_e,1) + 1.930000E+01;
        if(M_e >= 5.5){
          a1 = -9.703704E-03 * pow(M_e,3) +1.998571E-01* pow(M_e,2) -1.374931E+00* pow(M_e,1) + 3.180857E+00;
        }

        a2 = -9.603730E-03 * pow(M_e,4) +2.013313E-01* pow(M_e,3) -1.597797E+00* pow(M_e,2) + 5.974929E+00* pow(M_e,1) -9.764600E+00;

        a3 = 6.049036E-03 * pow(M_e,3) -5.273554E-02* pow(M_e,2) + 6.116542E-02* pow(M_e,1) + 2.744318E-01;
        if(M_e >= 5.5){
         a3 = -3.740741E-04 * pow(M_e,3) +1.293571E-02* pow(M_e,2) -1.296886E-01* pow(M_e,1) + 4.062429E-01;
        }

        a4 = 2.937939E-02 * pow(M_e,3) -4.928550E-01* pow(M_e,2) + 2.853076E+00* pow(M_e,1) -5.828876E+00;
        if(M_e >= 5.5){
         a4 = 4.586667E-02 * pow(M_e,3) -1.042600E+00* pow(M_e,2) + 7.882633E+00* pow(M_e,1) -1.982340E+01;
        }

        dNdRet = a1 * exp(a2 * H12) + a3 * exp(a4 * H12);
        */

        //ver 3.4
        /*
        a1 = -2.533333E-02 * pow(M_e,3) + 3.520000E-01 * pow(M_e, 2) - 1.677667E+00 * pow(M_e, 1) + 2.880000E+00;
        if(M_e >= 5.5){
          a1 = 1.333333E-03 * pow(M_e, 4) - 3.851852E-02 * pow(M_e,3) + 4.036667E-01 * pow(M_e, 2) - 1.823228E+00 * pow(M_e, 1) + 3.091286E+00;
        }

        a2 = 2.135237E-01 * pow(M_e,3) - 3.355008E+00 * pow(M_e, 2) + 1.756479E+01 * pow(M_e, 1) - 3.116920E+01;
        if(M_e >= 5.5){
         a2 = 3.963138E-02 * pow(M_e,3) - 8.250190E-01 * pow(M_e, 2) + 5.709339E+00 * pow(M_e, 1) - 1.356531E+01;
        }

        a3 = -4.117364E-03 * pow(M_e,3) + 6.827581E-02 * pow(M_e, 2) - 3.643231E-01 * pow(M_e, 1) + 6.363908E-01;
        if(M_e >= 5.5){
         a3 = 8.367309E-04 * pow(M_e,3) - 1.332571E-02 * pow(M_e, 2) + 6.164664E-02 * pow(M_e, 1) - 6.219784E-02;
        }

        a4 = 3.084084E-01 * pow(M_e,3) - 4.565485E+00 * pow(M_e, 2) + 2.246646E+01 * pow(M_e, 1) - 3.690620E+01;
        if(M_e >= 5.5){
         a4 = 2.938777E-03 * pow(M_e,3) - 9.799887E-02 * pow(M_e, 2) + 9.666538E-01 * pow(M_e, 1) - 2.976438E+00;
        }

        dNdRet = a1 * exp(a2 * H12) + a3 * exp(a4 * H12);
        //ver3.4.1
        su2double mindNdRet = -0.0003397773*pow(M_e,3) + 0.007076021*pow(M_e,2)- 0.05221073*pow(M_e,1) + 0.14173;
        */

        //ver3.4.2
        a1 = -4.000000E-02 *  pow(M_e,3) + 5.000000E-01 *  pow(M_e,2) - 1.900000E+00 *  pow(M_e,1) + 2.250000E+00;
        if(M_e >= 5.5){
          a1 = -1.833333E-02 *  pow(M_e,4) + 5.057407E-01 *  pow(M_e,3) - 5.187917E+00 *  pow(M_e,2) + 2.341089E+01 *  pow(M_e,1) - 3.892161E+01;
        }

        a2 = 3.333333E-02 *  pow(M_e,3) - 4.900000E-01 *  pow(M_e,2) + 2.416667E+00 *  pow(M_e,1) - 4.640000E+00;
        if(M_e >= 5.5){
         a2 = -3.333333E-03 *  pow(M_e,4) + 1.233333E-01 *  pow(M_e,3) - 1.549167E+00 *  pow(M_e,2) + 8.194167E+00 *  pow(M_e,1) - 1.630000E+01;
        }

        a3 = 8.000000E-04 *  pow(M_e,3) - 1.480000E-02 *  pow(M_e,2) + 9.040000E-02 *  pow(M_e,1) - 1.750000E-01;
        if(M_e >= 5.5){
         a3 = 1.500000E-04 *  pow(M_e,4) - 2.875926E-03 *  pow(M_e,3) + 1.768750E-02 *  pow(M_e,2) - 3.266852E-02 *  pow(M_e,1) - 6.562500E-03;
        }

        a4 = 2.666667E-02 *  pow(M_e,3) - 3.400000E-01 *  pow(M_e,2) + 1.463333E+00 *  pow(M_e,1) - 2.300000E+00;
        if(M_e >= 5.5){
         a4 = -5.500000E-03 *  pow(M_e,4) + 1.500926E-01 *  pow(M_e,3) - 1.544875E+00 *  pow(M_e,2) + 7.140423E+00 *  pow(M_e,1) - 1.257873E+01;
        }

        dNdRet = a1 * exp(a2 * H12) + a3 * exp(a4 * H12);
        //ver3.4.2
        su2double mindNdRet = 2.210826E-04*  pow(M_e,3) - 1.966442E-03*  pow(M_e,2) - 4.500020E-03*  pow(M_e,1) + 5.940342E-02;

        dNdRet = max(dNdRet, 0.00000000001);
        dNdRet = min(dNdRet, mindNdRet);
        break;
      }

      case AFMT_CORRELATION::DEFAULT:
        SU2_MPI::Error("Transition correlation is set to DEFAULT but no default value has ben set in the code.",
                       CURRENT_FUNCTION);
        break;
    }
    /*
    if(dNdRet == nan){
      cout << "NaN is tected." << endl;
    }
    */
    return dNdRet;
  }


  /*!
   * \brief Compute dN/dRet from correlations.
   * \param[in] H12 - Integreated Shape Factor.
   * \param[in] Hk - Kinetic Shape Factor.
   * \param[in] M_e - Edge Mach number.
   * \param[out] dNdRet - N factor gradient for Mack 2nd mode.
   */
  su2double Ret0_Correlations(const su2double H12, const su2double Hk, const su2double M_e) const {
    su2double Ret0 = 0.0;
    su2double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0;

    switch (options.Correlation) {
      case AFMT_CORRELATION::Liu2023: {
        a1 = 2.704e-1 *  pow(M_e,2) - 5.310 * M_e + 32.16;
        a2 = -1.423e-3 * pow(M_e,3) + 2.791e-2 * pow(M_e,2) - 1.610e-1 * M_e - 3.092;
        Ret0 = a1 * pow(Hk,a2) + 2.0;
        Ret0 = pow(10,Ret0);
        break;
      }

      case AFMT_CORRELATION::sok: {
        /*
        //ver1
        a1 = - 6.6660e-04 * pow(M_e,3) + 1.4219e-02 * pow(M_e,2) - 1.0373e-01 * M_e + 2.6474e-01;
        a2 = - 1.2833e-02 * pow(M_e,3) + 2.7246e-01 * pow(M_e,2) - 1.9261 * M_e + 6.3073;
        a3 = + 1.2615e-01 * pow(M_e,3) - 2.5264 * pow(M_e,2) + 1.7195e+01 * M_e - 3.9883e+01;
        */

        //ver2
        /*
        a1 = + 4498 * pow(M_e, - 7.966) + 0.007729;
        a2 = + 0.03268 * pow(M_e, 3) - 0.645 * pow(M_e, 2) + 4.071 * pow(M_e, 1) - 6.372;
        a3 = + 0.01727 * pow(M_e, 4) - 0.593 * pow(M_e, 3) + 7.099 * pow(M_e, 2) - 34.87 * M_e + 59.2;
        */

        //ver2.3 custom
        /*
        a1 = + 1.0320e-02 * pow(M_e,2) - 1.3351e-01 * M_e + 4.3280e-01;
        a2 = - 2.2190e-01 * pow(M_e,2) + 2.6323 * M_e - 5.6543;
        a3 = + 1.2667 * pow(M_e,2) - 1.4497e+01 * M_e + 3.8812e+01;
        */

        //ver3.0 custom
        a1 = - 7.8698E-04 * pow(M_e,3) + 1.8019E-02 * pow(M_e,2) - 1.3741E-01 * M_e + 3.5579E-01;
        a2 = + 8.8805E-03 * pow(M_e,3) - 1.7327E-01 * pow(M_e,2) + 1.0372E+00 * M_e - 6.0120E-02;
        a3 = + 8.2566E-03 * pow(M_e,3) - 1.4756E-01 * pow(M_e,2) + 1.6266E+00 * M_e - 6.5092E+00;


        Ret0 = (a1 * pow(H12,2) + a2 * H12 + a3) / (H12 - 1.5);
        Ret0 = max(1.9, min(Ret0,5.0));
        Ret0 = pow(10, Ret0);
        break;
      }

      case AFMT_CORRELATION::DEFAULT:
        SU2_MPI::Error("Transition correlation is set to DEFAULT but no default value has ben set in the code.",
                       CURRENT_FUNCTION);
        break;
    }
    /*
    if(Ret0 == nan){
      cout << "NaN is tected." << endl;
    }
    */
    return Ret0;
  }




  /*!
   * \brief Compute D_H12_y/theta from correlations.
   * \param[in] H12 - Integreated Shape Factor.   
   * \param[in] M_e - Edge Mach number.
   * \param[out] dNdRet - N factor gradient for Mack 2nd mode.
   */
  su2double D_H12_Correlations(const su2double H12, const su2double Hk) const {
    su2double D_H12 = 0.0;
    su2double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, a5 = 0.0, a6 = 0.0, a7 =0.0;

    switch (options.Correlation) {
      case AFMT_CORRELATION::Liu2023: {
        /*a1 = 9.1610e-4;
        a2 = -2.9340e-3;
        a3 = 1.6520;
        a4 = 8.9230e-1;
        a5 = -6.6490;
        a6 = -3.6160e-1;
        a7 = 1.3950e+1; 
        D_H12 = a1 * pow(H12, 2) * Hk + a2 * pow(H12, 2) + a3 * pow(H12, 1) + a4 * pow(Hk, 2) + a5 * pow(Hk, 1) + a6 * H12 * Hk + a7;
        */
        D_H12 = 0.4572 * H12 + 0.386 * Hk + 0.04384 * H12 * Hk - 0.0002832 * pow(H12,2) + 1.164;
        break;
      }

      case AFMT_CORRELATION::sok: {
        D_H12 = 0.4572 * H12 + 0.386 * Hk + 0.04384 * H12 * Hk - 0.0002832 * pow(H12,2) + 1.164;
        break;
      }

      case AFMT_CORRELATION::DEFAULT:
        SU2_MPI::Error("Transition correlation is set to DEFAULT but no default value has ben set in the code.",
                       CURRENT_FUNCTION);
        break;
    }
    D_H12 = min(40.0, max(1.0, D_H12));
    return D_H12;
  }


  /*!
   * \brief Compute l_H12_y/theta from correlations.
   * \param[in] H12 - Integreated Shape Factor.   
   * \param[in] Hk - Edge Mach number.
   * \param[out] l(H12,Hk) - l function.
   */
  su2double l_H12_Correlations(const su2double H12, const su2double Hk) const {
  su2double l_H12 = 0.0;
  su2double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, a5 = 0.0, a6 = 0.0, a7 =0.0;

    switch (options.Correlation) {
      case AFMT_CORRELATION::Liu2023: {
        /*
        a1 = 9.1610e-4;
        a2 = -2.9340e-3;
        a3 = 1.6520;
        a4 = 8.9230e-1;
        a5 = -6.6490;
        a6 = -3.6160e-1;
        a7 = 1.3950e+1; 
        l_H12 = a1 * pow(H12, 2) * Hk + a2 * pow(H12, 2) + a3 * pow(H12, 1) + a4 * pow(Hk, 2) + a5 * pow(Hk, 1) + a6 * H12 * Hk + a7;
        */
        l_H12 = 0.1529 - 0.002641 * H12 + 0.2895 * Hk + 0.0005796 * pow(H12,2) - 0.01232 * H12 * Hk - 0.06548 * pow(Hk,2) - 8.154e-07 * pow(H12,3) - 0.0001493 * pow(H12,2) * Hk + 0.003404 * H12 * pow(Hk,2); //polyfitting
        break;
      }

      case AFMT_CORRELATION::sok: {
        l_H12 = 0.1529 - 0.002641 * H12 + 0.2895 * Hk + 0.0005796 * pow(H12,2) - 0.01232 * H12 * Hk - 0.06548 * pow(Hk,2) - 8.154e-07 * pow(H12,3) - 0.0001493 * pow(H12,2) * Hk + 0.003404 * H12 * pow(Hk,2); //polyfitting
        break;
      }

      case AFMT_CORRELATION::DEFAULT:
        SU2_MPI::Error("Transition correlation is set to DEFAULT but no default value has ben set in the code.",
                       CURRENT_FUNCTION);
        break;
    }
    l_H12 = min(0.6, max(0.1, l_H12));
    return l_H12;
  }



};
