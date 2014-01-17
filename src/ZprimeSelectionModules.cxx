#include "include/ZprimeSelectionModules.h"

IsoConeSelection::IsoConeSelection(TString type, double a, double b, double c, double iso_num){
 m_type = type;
 m_a = a;
 m_b = b;
 m_c = c;
 m_iso = iso_num;

}

bool IsoConeSelection::pass(BaseCycleContainer *bcc){

 double pt;
 Particle particle;

 //needs adaptation for the electron channel 
 if(m_type.Contains("mu")||m_type.Contains("Mu")){  
   pt = bcc->muons->at(0).pt();
   particle = bcc->muons->at(0);
 }
 else{
   return false;
 }

 double radius = m_a/(pt+m_b)+m_c; 
 if(radius <0.02) radius = 0.02;
 double value = relIso(particle,radius); 


 if(relIso(particle,0.4) < 0.2) return true;

 return  value < m_iso? true : false;
}

std::string IsoConeSelection::description()
{
    char s[100];
    sprintf(s, "IsoCone parameters:  a = %f, b = %f, c = %f; Iso = %f",m_a,m_b,m_c,m_iso);
    return s;
}





