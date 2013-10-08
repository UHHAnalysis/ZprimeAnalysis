// Dear emacs, this is -*- c++ -*-
#ifndef TopFitCalc_H
#define TopFitCalc_H

#include "EventCalc.h"
#include "TVector3.h"
#include "Utils.h"
#include "TTbarGen.h"
#include "TMinuit.h"

/**
 * class for the calculation of Tops with hadronic Tags
 *
 * author Daniel Gonzalez
 */

class TopFitCalc
{
 private:
  static TopFitCalc* m_instance;

 public:
  static TopFitCalc* Instance();

  /*
   * Take the C/A-Jet with the highest pT and use it as hadronic Top. 
   * For the leptonic side exactly one muon is expected. A list is generated 
   * with all possible combinations of one ak 0.5 Jets. Later one ak 0.5 Jet is selected  
   * via chi_square.
   * ak 0.5 Jets with  
   */
  void CalculateSelection();

  /*
   * Take the C/A-Jet with the TopTag and use it as hadronic Top.
   * For the leptonic side exactly one muon is expected. A list is generated 
   * with all possible combinations of one ak 0.5 Jets. Later one ak 0.5 Jet is selected 
   * via chi_square. Make sure a reasonable selection is run. There is hardcoded a minimum
   * distance of muon and TopTag of 0.8 and between ak 0.5 and TopTag of 1.3 in delta R
   *
   */
  void CalculateTopTag();


  void Reset();


  double delPhi(double x1, double x2) const{// x1 & x2 are two Phi expected in the range [-PI,PI]
    double deltaphi = fabs(x1 - x2);
    if(deltaphi > M_PI) deltaphi = 2*M_PI - deltaphi;
    return deltaphi;
  }

  double delR(LorentzVector vector1, LorentzVector vector2){
   double phi = delPhi(vector1.phi(),vector2.phi());
   double eta = vector1.eta()-vector2.eta();	
   return sqrt(phi*phi+eta*eta);
  }	



  /*
   * Basic idea is that for the case the Neutrino becomes an Im-Part you do a fit
   * on the px and py components. If you look into the kinematics there is a way of 
   * doing this with a W-mass constraint. As the name indicates polarcoordinates 
   * are used in the computation.
   *
   */
  std::vector<LorentzVector> NeutrinoFitPolar(const LorentzVector lepton, const LorentzVector met);

 private:

  mutable SLogger m_logger;
  TopFitCalc();
  ~TopFitCalc();


  TMinuit* positiv;// For the Neutrino-fit

  // data members to store calculated results

  bool b_Reconstruction;

};



#endif // TopFitCalc
