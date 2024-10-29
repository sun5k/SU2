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
        /*
        H12 = (1.0696e+1 * pow(T_over_T0,2) + 1.4386e+1 * T_over_T0 * M_e - 2.7767e+2 * T_over_T0 - 1.8548e+1 * M_e + 2.2921e+2) * T_over_T0;
        H12 = H12 + 1.7856 * pow(M_e,2) + 5.9571 * M_e - 5.7666e+1;
        */
        H12 = (1.7224E+02 * pow(T_over_T0, 3) + (2.0322E+01 * pow(T_over_T0, 2) * M_e)) + (-4.1141E+02 * pow(T_over_T0, 2));
        H12 = H12 + (1.7533E+00 * pow(M_e, 2)) + (-2.4781E+01 * T_over_T0 * M_e);
        H12 = H12 + (3.1445E+02 * T_over_T0) + (7.8102E+00 * M_e) - 7.5219E+01;
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
        /*
        a1 = 2.642158e-5 * pow(M_e,5) - 1.047983e-3 * pow(M_e,4) + 1.688650e-2 * pow(M_e,3) - 1.392065e-1 * pow(M_e,2) + 5.936119e-1 * M_e - 1.065709;
        a2 = 1.140608e-4 * pow(M_e,4) - 9.542641e-3 * pow(M_e,3) + 2.129670e-1 * pow(M_e,2) - 1.968956 * M_e + 7.184406;
        a3 = 4.129277e-3 * pow(M_e,4) - 1.256703e-1 * pow(M_e,3) + 1.403225 * pow(M_e,2) - 6.410845 * M_e + 1.309789e+1;
        */
        a1 = -6.3666E-01 * exp(-7.1793E-01 * M_e) + -5.2614E-02 * exp(-3.0234E-01 * M_e);
        a2 = +2.9624E-02 * pow(M_e, 2) - 6.0021E-01 * M_e + 3.5517E+00;
        a3 = 8.6991E+00 * exp(-4.6608E-01 * M_e) + 3.4913E+00 * exp(8.0897E-03 * M_e);
        Hk = a1 * pow(H12,2) + a2 * H12 + a3;
        Hk = log(Hk);
        break;
      }

      case AFMT_CORRELATION::sok: {
        /*
        a1 = 2.642158e-5 * pow(M_e,5) - 1.047983e-3 * pow(M_e,4) + 1.688650e-2 * pow(M_e,3) - 1.392065e-1 * pow(M_e,2) + 5.936119e-1 * M_e - 1.065709;
        a2 = 1.140608e-4 * pow(M_e,4) - 9.542641e-3 * pow(M_e,3) + 2.129670e-1 * pow(M_e,2) - 1.968956 * M_e + 7.184406;
        a3 = 4.129277e-3 * pow(M_e,4) - 1.256703e-1 * pow(M_e,3) + 1.403225 * pow(M_e,2) - 6.410845 * M_e + 1.309789e+1;
        */
        a1 = -6.3666E-01 * exp(-7.1793E-01 * M_e) + -5.2614E-02 * exp(-3.0234E-01 * M_e);
        a2 = +2.9624E-02 * pow(M_e, 2) - 6.0021E-01 * M_e + 3.5517E+00;
        a3 = 8.6991E+00 * exp(-4.6608E-01 * M_e) + 3.4913E+00 * exp(8.0897E-03 * M_e);
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
  su2double RevRet_Correlations(const su2double H12, const su2double M_e, const su2double T_e) const {
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
        /*
        a1 = 0.0008 * pow(M_e,2) + 0.091 * M_e + 0.1253;
        a2 = -0.0343 * pow(M_e,2) - 0.1322 * M_e + 0.8943;
        a3 = 0.0887 * pow(M_e,2) - 0.6799 * M_e + 2.3288;
        RevRet = a1 * pow(log(H12),2) + a2 * log(H12) + a3;
        */
        RevRet = -1.4461e-03 * pow(H12, 2) + 1.1681e-07 * H12 * pow(T_e, 2) - 9.8006e-05 * H12 * T_e + 2.3730e-01 * H12 ;
        RevRet = RevRet - 3.5432e-01 * M_e- 1.9925e-06 * pow(T_e, 2) + 1.8258e-03 * T_e + 2.2939e+00;
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
    RevRet = max(min(10.0, RevRet),1.0);
    return RevRet;
  }





  /*!
   * \brief Compute dN/dRet from correlations.
   * \param[in] H12 - Integreated Shape Factor.   
   * \param[in] M_e - Edge Mach number.
   * \param[out] dNdRet - N factor gradient for Mack 2nd mode.
   */
  su2double dNdRet_Correlations(const su2double H12, const su2double M_e, const su2double T_e) const {
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
        
       //ver 3.6
        /*
        a1 = 1.333333E-02 * pow(M_e, 3) - 3.400000E-01 * pow(M_e, 2) + 2.426667E+00 * pow(M_e, 1) - 5.070000E+00;
        if(M_e >= 5.5){
          a1 = -5.953298E-03 * pow(M_e, 4) + 1.776048E-01 * pow(M_e, 3) - 1.939709E+00 * pow(M_e, 2) + 9.124382E+00 * pow(M_e, 1) - 1.539809E+01;
        }

        a2 = 5.333333E-02 * pow(M_e, 3) - 7.400000E-01 * pow(M_e, 2) + 3.536667E+00 * pow(M_e, 1) - 6.620000E+00;
        if(M_e >= 5.5){
         a2 = -2.475387E-02 * pow(M_e, 4) + 6.806654E-01 * pow(M_e, 3) - 7.027497E+00 * pow(M_e, 2) + 3.240888E+01 * pow(M_e, 1) - 5.694479E+01;
        }

        a3 = 9.333333E-03 * pow(M_e, 3) - 1.220000E-01 * pow(M_e, 2) + 4.986667E-01 * pow(M_e, 1) - 5.900000E-01;
        if(M_e >= 5.5){
         a3 = 1.299084E-03 * pow(M_e, 4) - 3.829186E-02 * pow(M_e, 3) + 4.215650E-01 * pow(M_e, 2) - 2.054308E+00 * pow(M_e, 1) + 3.743436E+00;
        }

        a4 = 1.327229E-01 * pow(M_e, 3) - 2.011760E+00 * pow(M_e, 2) + 1.035974E+01 * pow(M_e, 1) - 1.834507E+01;
        if(M_e >= 5.5){
         a4 = -1.038156E-02 * pow(M_e, 4) + 3.088321E-01 * pow(M_e, 3) - 3.438725E+00 * pow(M_e, 2) + 1.698829E+01 * pow(M_e, 1) - 3.143728E+01;
        }


        dNdRet = a1 * exp(a2 * H12) + a3 * exp(a4 * H12);
        su2double mindNdRet = 8.411474E-05 * pow(M_e, 3) - 9.546653E-04 * pow(M_e, 2) - 1.609030E-03 * pow(M_e, 1) + 3.550163E-02;
        */

       //ver 4.0
       /*
        a1 = 1.4141e-02 * pow(M_e, 3) + 1.0476e-04 * pow(M_e, 2) * T_e - 2.5729e-01 * pow(M_e, 2) - 1.3143e-03 * M_e * T_e + 1.4308e+00 * M_e + 4.1345e-03 * T_e - 2.1285e+00;
        a2 = 1.0909e-02 * pow(M_e, 3) - 1.0952e-04 * pow(M_e, 2) * T_e - 2.7145e-01 * pow(M_e, 2) + 1.3833e-03 * M_e * T_e + 2.2240e+00 * M_e - 4.0988e-03 * T_e - 6.6034e+00;
        a3 = 6.1515e-04 * pow(M_e, 4) - 4.0404e-07 * pow(M_e, 3) * T_e - 1.7920e-02 * pow(M_e, 3) + 6.7186e-06 * pow(M_e, 2) * T_e + 1.9604e-01 * pow(M_e, 2);
        a3 = a3 - 3.6700e-05 * M_e * T_e - 9.5263e-01 * M_e + 6.5946e-05 * T_e + 1.7377e+00;
        a4 = -6.9077e-03 * pow(M_e, 5) + 1.9697e-06 * pow(M_e, 4) * T_e + 2.1903e-01 * pow(M_e, 4) - 4.8232e-05 * pow(M_e, 3) * T_e - 2.7339e+00 * pow(M_e, 3);
        a4 = a4 + 4.3617e-04 * pow(M_e, 2) * T_e + 1.6721e+01 * pow(M_e, 2) - 1.7270e-03 * M_e * T_e - 4.9803e+01 * M_e + 2.5271e-03 * T_e + 5.7184e+01;

        if ( T_e < 200){
          a1 = 1.1212e-02 * pow(M_e, 3) - 1.6667e-04 * pow(M_e, 2) * T_e - 1.4808e-01 * pow(M_e, 2) + 2.4278e-03 * M_e * T_e + 3.4589e-01 * M_e - 9.0861e-03 * T_e + 1.1884e+00;
          a2 = 1.9394e-02 * pow(M_e, 3) + 4.2539e-04 * pow(M_e, 2) * T_e - 5.3753e-01 * pow(M_e, 2) - 5.7429e-03 * M_e * T_e + 4.6240e+00 * M_e + 1.8867e-02 * T_e - 1.3145e+01;
          a3 = 2.1894e-04 * pow(M_e, 4) - 3.0020e-05 * pow(M_e, 3) * T_e - 2.0917e-03 * pow(M_e, 3) + 6.0709e-04 * pow(M_e, 2) * T_e - 1.5627e-02 * pow(M_e, 2) ;
          a3 = a3 - 4.0613e-03 * M_e * T_e + 2.2338e-01 * M_e + 8.9619e-03 * T_e - 5.9703e-01;
          a4 = -3.3154e-03 * pow(M_e, 5) + 7.1566e-05 * pow(M_e, 4) * T_e + 9.3194e-02 * pow(M_e, 4) - 1.7594e-03 * pow(M_e, 3) * T_e - 1.0113e+00 * pow(M_e, 3) ;
          a4 = a4 + 1.5972e-02 * pow(M_e, 2) * T_e + 5.1928e+00 * pow(M_e, 2) - 6.3218e-02 * M_e * T_e - 1.2095e+01 * M_e + 9.1490e-02 * T_e + 9.0577e+00;
        }
        */


       //ver 4.1
       /*
       a1 = -5.0370e-04 * pow(M_e, 3) * T_e + 1.2519e-01 * pow(M_e, 3) + 1.0386e-02 * pow(M_e, 2) * T_e - 2.5071e+00 * pow(M_e, 2) - 7.0567e-02 * M_e * T_e + 1.6451e+01 * M_e + 1.5792e-01 * T_e - 3.5117e+01;
       a2 = 4.0741e-04 * pow(M_e, 3) * T_e - 8.7407e-02 * pow(M_e, 3) - 8.1571e-03 * pow(M_e, 2) * T_e + 1.6514e+00 * pow(M_e, 2) + 5.3605e-02 * M_e * T_e - 1.0092e+01 * M_e - 1.1516e-01 * T_e + 1.9136e+01;
       a3 = -7.3333e-06 * pow(M_e, 4) * T_e + 2.9333e-03 * pow(M_e, 4) + 1.9622e-04 * pow(M_e, 3) * T_e - 8.0370e-02 * pow(M_e, 3) - 1.9572e-03 * pow(M_e, 2) * T_e + 8.2297e-01 * pow(M_e, 2) + 8.6217e-03 * M_e * T_e - 3.7325e+00 * M_e - 1.4151e-02 * T_e + 6.3334e+00;
       a4 = 2.3333e-05 * pow(M_e, 4) * T_e - 1.8333e-02 * pow(M_e, 4) - 6.3519e-04 * pow(M_e, 3) * T_e + 5.1226e-01 * pow(M_e, 3) + 6.4792e-03 * pow(M_e, 2) * T_e - 5.3659e+00 * pow(M_e, 2) - 2.9326e-02 * M_e * T_e + 2.4978e+01 * M_e + 4.9623e-02 * T_e - 4.3609e+01;
       
       if (M_e < 5.5) {
        a1 = 2.0000e-04 * pow(M_e, 2) * T_e - 6.0000e-02 * pow(M_e, 2) - 2.3000e-03 * M_e * T_e + 5.9000e-01 * M_e + 6.8000e-03 * T_e - 1.0900e+00;
        a2 = 2.4000e-03 * pow(M_e, 2) * T_e - 3.2000e-01 * pow(M_e, 2) - 2.2800e-02 * M_e * T_e + 3.2400e+00 * M_e + 5.3500e-02 * T_e - 9.1000e+00;
        a3 = 1.0000e-04 * pow(M_e, 2) * T_e - 1.6000e-02 * pow(M_e, 2) - 1.1500e-03 * M_e * T_e + 1.6200e-01 * M_e + 3.3000e-03 * T_e - 3.9500e-01;
        a4 = 6.0000e-04 * pow(M_e, 2) * T_e - 2.2000e-01 * pow(M_e, 2) - 6.0000e-03 * M_e * T_e + 2.4800e+00 * M_e + 1.4850e-02 * T_e - 7.0850e+00;
       }

       if ( T_e < 200 ) {
        a1 = 1.9753e-04 * pow(M_e, 3) * T_e - 1.5062e-02 * pow(M_e, 3) - 3.8286e-03 * pow(M_e, 2) * T_e + 3.3571e-01 * pow(M_e, 2) + 2.4589e-02 * M_e * T_e - 2.5803e+00 * M_e - 5.2743e-02 * T_e + 7.0157e+00;
        a2 = -2.6667e-04 * pow(M_e, 3) * T_e + 4.7407e-02 * pow(M_e, 3) + 5.0190e-03 * pow(M_e, 2) * T_e - 9.8381e-01 * pow(M_e, 2) - 3.1081e-02 * M_e * T_e + 6.8448e+00 * M_e + 6.1995e-02 * T_e - 1.6295e+01;
        a3 = 1.2000e-05 * pow(M_e, 4) * T_e - 9.3333e-04 * pow(M_e, 4) - 3.2775e-04 * pow(M_e, 3) * T_e + 2.4425e-02 * pow(M_e, 3) + 3.3930e-03 * pow(M_e, 2) * T_e - 2.4707e-01 * pow(M_e, 2) - 1.5751e-02 * M_e * T_e + 1.1421e+00 * M_e + 2.7523e-02 * T_e - 2.0013e+00;
        a4 = -7.5556e-05 * pow(M_e, 4) * T_e + 1.4444e-03 * pow(M_e, 4) + 2.1477e-03 * pow(M_e, 3) * T_e - 4.4309e-02 * pow(M_e, 3) - 2.2881e-02 * pow(M_e, 2) * T_e + 5.0614e-01 * pow(M_e, 2) + 1.0828e-01 * M_e * T_e - 2.5442e+00 * M_e - 1.9198e-01 * T_e + 4.7124e+00;

        if (M_e < 5.5) {
          a1 = 5.3333e-04 * pow(M_e, 2) * T_e - 1.2667e-01 * pow(M_e, 2) - 4.8000e-03 * M_e * T_e + 1.0900e+00 * M_e + 9.8000e-03 * T_e - 1.6900e+00;
          a2 = -2.6667e-04 * pow(M_e, 2) * T_e + 2.1333e-01 * pow(M_e, 2) - 8.0000e-04 * M_e * T_e - 1.1600e+00 * M_e + 1.1000e-02 * T_e - 6.0000e-01;
          a3 = -5.3333e-05 * pow(M_e, 2) * T_e + 1.4667e-02 * pow(M_e, 2) + 2.6667e-04 * M_e * T_e - 1.2133e-01 * M_e + 1.2667e-04 * T_e + 2.3967e-01;
          a4 = 1.0667e-03 * pow(M_e, 2) * T_e - 3.1333e-01 * pow(M_e, 2) - 1.0200e-02 * M_e * T_e + 3.3200e+00 * M_e + 2.3433e-02 * T_e - 8.8017e+00;
        }
       }
       */

      //ver 4.2

       a1 = 1.7037e-04 * pow(M_e, 3) * T_e + 3.7037e-03 * pow(M_e, 3) - 3.6857e-03 * pow(M_e, 2) * T_e + 4.5714e-02 * pow(M_e, 2) + 2.6186e-02 * M_e * T_e - 1.2481e+00 * M_e - 6.0850e-02 * T_e + 5.3279e+00;
       a2 = -2.8519e-04 * pow(M_e, 3) * T_e + 5.5556e-02 * pow(M_e, 3) + 5.4000e-03 * pow(M_e, 2) * T_e - 1.0886e+00 * pow(M_e, 2) - 3.3593e-02 * M_e * T_e + 7.1547e+00 * M_e + 6.8946e-02 * T_e - 1.6464e+01;
       a3 = 1.2833e-05 * pow(M_e, 4) * T_e - 2.5833e-03 * pow(M_e, 4) - 3.4335e-04 * pow(M_e, 3) * T_e + 6.7872e-02 * pow(M_e, 3) + 3.4200e-03 * pow(M_e, 2) * T_e - 6.6143e-01 * pow(M_e, 2) - 1.5029e-02 * M_e * T_e + 2.8305e+00 * M_e + 2.4587e-02 * T_e - 4.4778e+00;
       a4 = 2.1333e-05 * pow(M_e, 5) * T_e - 3.4667e-03 * pow(M_e, 5) - 7.7867e-04 * pow(M_e, 4) * T_e + 1.2093e-01 * pow(M_e, 4) + 1.1240e-02 * pow(M_e, 3) * T_e - 1.6563e+00 * pow(M_e, 3) - 8.0253e-02 * pow(M_e, 2) * T_e + 1.1108e+01 * pow(M_e, 2) + 2.8355e-01 * M_e * T_e - 3.6346e+01 * M_e - 3.9675e-01 * T_e + 4.6132e+01;

       if (M_e < 5.5) {
        a1 = 2.0000e-04 * pow(M_e, 2) * T_e - 2.2000e-01 * pow(M_e, 2) - 2.5000e-03 * M_e * T_e + 2.1700e+00 * M_e + 7.7000e-03 * T_e - 4.8100e+00;
        a2 = 2.0000e-04 * pow(M_e, 2) * T_e - 1.2000e-01 * pow(M_e, 2) - 1.7000e-03 * M_e * T_e + 1.5600e+00 * M_e + 3.3000e-03 * T_e - 5.7300e+00;
        a3 = 4.0000e-05 * pow(M_e, 2) * T_e + 1.2000e-02 * pow(M_e, 2) - 4.2000e-04 * M_e * T_e - 1.4600e-01 * M_e + 1.1000e-03 * T_e + 4.5000e-01;
        a4 = -9.0000e-02 * pow(M_e, 2) + 1.1750e+00 * M_e - 3.8450e+00;
       }

       if ( T_e < 200 ) {
        a1 = 1.7778e-04 * pow(M_e, 3) * T_e + 2.2222e-03 * pow(M_e, 3) - 4.6571e-03 * pow(M_e, 2) * T_e + 2.4000e-01 * pow(M_e, 2) + 3.7875e-02 * M_e * T_e - 3.5858e+00 * M_e - 9.8567e-02 * T_e + 1.2871e+01;
        a2 = 4.4444e-05 * pow(M_e, 3) * T_e - 1.0370e-02 * pow(M_e, 3) - 8.5714e-04 * pow(M_e, 2) * T_e + 1.6286e-01 * pow(M_e, 2) + 5.1889e-03 * M_e * T_e - 6.0169e-01 * M_e - 9.9381e-03 * T_e - 6.8667e-01;
        a3 = -7.8889e-06 * pow(M_e, 4) * T_e + 1.5611e-03 * pow(M_e, 4) + 2.0888e-04 * pow(M_e, 3) * T_e - 4.2573e-02 * pow(M_e, 3) - 2.0292e-03 * pow(M_e, 2) * T_e + 4.2842e-01 * pow(M_e, 2) + 8.5593e-03 * M_e * T_e - 1.8872e+00 * M_e - 1.3240e-02 * T_e + 3.0875e+00;
        a4 = -2.1556e-05 * pow(M_e, 4) * T_e - 3.4889e-03 * pow(M_e, 4) + 7.2259e-04 * pow(M_e, 3) * T_e + 8.4259e-02 * pow(M_e, 3) - 8.9203e-03 * pow(M_e, 2) * T_e - 7.3069e-01 * pow(M_e, 2) + 4.8179e-02 * M_e * T_e + 2.6442e+00 * M_e - 9.6213e-02 * T_e - 3.2561e+00;


        if (M_e < 5.5) {
          a1 = -1.2000e-03 * pow(M_e, 2) * T_e + 6.0000e-02 * pow(M_e, 2) + 1.2467e-02 * M_e * T_e - 8.2333e-01 * M_e - 3.3800e-02 * T_e + 3.4900e+00;
          a2 = 1.8948e-18 * pow(M_e, 2) * T_e - 8.0000e-02 * pow(M_e, 2) - 2.5333e-03 * M_e * T_e + 1.7267e+00 * M_e + 1.4067e-02 * T_e - 7.8833e+00;
          a3 = 7.6000e-05 * pow(M_e, 2) * T_e + 4.8000e-03 * pow(M_e, 2) - 1.0207e-03 * M_e * T_e - 2.5867e-02 * M_e + 3.3013e-03 * T_e + 9.7333e-03;
          a4 = 1.3067e-03 * pow(M_e, 2) * T_e - 3.5133e-01 * pow(M_e, 2) - 1.2507e-02 * M_e * T_e + 3.6763e+00 * M_e + 2.8687e-02 * T_e - 9.5823e+00;

        }
       }








        dNdRet = a1 * exp(a2 * H12) + a3 * exp(a4 * H12);
        su2double mindNdRet = -3.3027e-04 * pow(M_e, 3) + 1.1678e-05 * pow(M_e, 2) * T_e + 4.4909e-03 * pow(M_e, 2) - 2.2004e-07 * M_e * pow(T_e, 2) ;
        mindNdRet = mindNdRet - 4.6365e-05 * M_e * T_e - 3.0232e-02 * M_e + 1.8735e-06 * pow(T_e, 2) - 4.2550e-04 * T_e + 1.3489e-01;

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
  su2double D_H12_Correlations(const su2double H12, const su2double Hk, const su2double T_e, const su2double M_e) const {
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
        //D_H12 = 0.4572 * H12 + 0.386 * Hk + 0.04384 * H12 * Hk - 0.0002832 * pow(H12,2) + 1.164;
        //D_H12 = 4.8995E-05 * pow(H12, 2) + -1.4000E-03 * H12 * Hk + 5.9840E-01 * H12 + 5.6240E-01 * Hk + 3.8493E-01;
        D_H12 = -6.0221e-05 * H12 * T_e + 6.0985e-01 * H12 - 9.9217e-02 * M_e - 3.9482e-04 * T_e + 2.6637e+00;
        break;
      }

      case AFMT_CORRELATION::DEFAULT:
        SU2_MPI::Error("Transition correlation is set to DEFAULT but no default value has ben set in the code.",
                       CURRENT_FUNCTION);
        break;
    }
    D_H12 = min(40.0, max(3.0, D_H12));
    return D_H12;
  }


  /*!
   * \brief Compute l_H12_y/theta from correlations.
   * \param[in] H12 - Integreated Shape Factor.   
   * \param[in] Hk - Edge Mach number.
   * \param[in] Te - Edge Temperature.
   * \param[out] l(H12,Hk) - l function.
   */
  su2double l_H12_Correlations(const su2double H12, const su2double Hk, const su2double T_e) const {
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
        /*
        l_H12 = 0.1529 - 0.002641 * H12 + 0.2895 * Hk + 0.0005796 * pow(H12,2) - 0.01232 * H12 * Hk - 0.06548 * pow(Hk,2) - 8.154e-07 * pow(H12,3) - 0.0001493 * pow(H12,2) * Hk + 0.003404 * H12 * pow(Hk,2); //polyfitting
        l_H12 =-2.4791E-04 * pow(H12,2)*Hk + 8.0177E-04 *pow(H12,2) + 6.7812E-06 * H12*Te*Hk + 7.9735E-03 * H12*Hk + -1.9448E-05 *H12*Te + -3.1423E-02 * H12;
        l_H12 = l_H12 + -1.4847E-04 * Hk*Te + -5.5737E-02 * Hk + - 5.3318E-10 * pow(Te,3) + 9.2202E-07 * pow(Te,2) + -1.5883E-04 * Te + 6.7869E-01;
        */
        if (T_e < 180) {
          l_H12 = 4.4574e-05 * pow(H12, 2) - 8.3713e-06 * H12 * T_e - 5.1466e-03 * H12 - 3.5782e-04 * Hk * T_e + 3.6538e-02 * Hk + 9.9770e-06 * pow(T_e, 2) - 2.0356e-03 * T_e + 5.6516e-01;
        }
        else {
          l_H12 = -2.4787e-04 * pow(H12, 2) * Hk + 8.0165e-04 * pow(H12, 2) + 6.7680e-06 * H12 * Hk * T_e + 7.9773e-03 * H12 * Hk - 1.9412e-05 * H12 * T_e; 
          l_H12 = l_H12 - 3.1433e-02 * H12 - 1.4818e-04 * Hk * T_e - 5.5839e-02 * Hk + 2.4221e-07 * pow(T_e, 2) + 1.0472e-04 * T_e + 6.4849e-01;
        }
        l_H12 = min(0.5, max(0.1, l_H12));
        break;
      }

      case AFMT_CORRELATION::DEFAULT:
        SU2_MPI::Error("Transition correlation is set to DEFAULT but no default value has ben set in the code.",
                       CURRENT_FUNCTION);
        break;
    }
    l_H12 = min(0.5, max(0.1, l_H12));
    return l_H12;
  }



};
