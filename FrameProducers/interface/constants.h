#ifndef RecoE2E_constants_h
#define RecoE2E_constants_h

static const unsigned int Nproj = 5;
static const unsigned int Nadjproj = 2;
static const unsigned int Nhitproj = 2;
static const unsigned int nDetFrameH = 280;
static const unsigned int nDetFrameW = 360;

static const bool debug = false;

static const std::string adj_projections[Nadjproj] = {"_5x5", "_3x3"};
static const int HBHE_IPHI_NUM = hcaldqm::constants::IPHI_NUM;//72;
unsigned int granularityMultiPhi[Nadjproj];
unsigned int granularityMultiEta[Nadjproj];
      
int totalEtaBins[Nadjproj];// = totalMultiEta*(eta_nbins_HBHE);
int totalPhiBins[Nadjproj];// = granularityMultiPhi * granularityMultiECAL*HBHE_IPHI_NUM;
std::vector<double> adjEtaBins[Nadjproj];

static const int nTOB = 6;
static const int nTEC = 9;
static const int nTIB = 4;
static const int nTID = 3;
static const int nBPIX = 4;
static const int nFPIX = 3;
static const int nEE = 2;
static const float zs = 0.;

static const int EB_IPHI_MIN = EBDetId::MIN_IPHI;//1;
static const int EB_IPHI_MAX = EBDetId::MAX_IPHI;//360;
static const int EB_IETA_MIN = EBDetId::MIN_IETA;//1;
static const int EB_IETA_MAX = EBDetId::MAX_IETA;//85;
static const int EE_MIN_IX = EEDetId::IX_MIN;//1;
static const int EE_MIN_IY = EEDetId::IY_MIN;//1;
static const int EE_MAX_IX = EEDetId::IX_MAX;//100;
static const int EE_MAX_IY = EEDetId::IY_MAX;//100;
static const int EE_NC_PER_ZSIDE = EEDetId::IX_MAX*EEDetId::IY_MAX; // 100*100
static const int HBHE_IETA_MAX_FINE = 20;
static const int HBHE_IETA_MAX_HB = hcaldqm::constants::IETA_MAX_HB;//16;
static const int HBHE_IETA_MIN_HB = hcaldqm::constants::IETA_MIN_HB;//1;
static const int HBHE_IETA_MAX_HE = hcaldqm::constants::IETA_MAX_HE;//29;
static const int HBHE_IETA_MAX_EB = hcaldqm::constants::IETA_MAX_HB + 1; // 17
static const int HBHE_IPHI_NUM = hcaldqm::constants::IPHI_NUM;//72;
static const int HBHE_IPHI_MIN = hcaldqm::constants::IPHI_MIN;//1;
static const int HBHE_IPHI_MAX = hcaldqm::constants::IPHI_MAX;//72;
static const int ECAL_IETA_MAX_EXT = 140;

static const std::string projections[Nproj] = {"", "_atECAL", "_atHCAL","_atECALfixIP","_atECALfixIPfromPV"}; //57425
static const std::string hit_projections[Nhitproj] = {"", "_atPV"};
static const int eta_nbins_HBHE = 2*(HBHE_IETA_MAX_HE-1);
static const int granularityMultiECAL=5;

// HBHE eta bin edges
static const double eta_bins_HBHE[2*(hcaldqm::constants::IETA_MAX_HE-1)+1] =
                  {-3.000, -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305,
                   -1.218, -1.131, -1.044, -0.957, -0.870, -0.783, -0.695, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0.000,
                    0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,  0.695,  0.783,  0.870,  0.957,  1.044,  1.131,  1.218,
                    1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,  1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  3.000}; // 57
//
// static data member definitions
//

// EE-(phi,eta) projection eta edges
// These are generated by requiring 5 fictional crystals
// to uniformly span each HCAL tower in eta (as in EB).
static const double eta_bins_EEm[5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB)+1] =
                  {-3.    , -2.93  , -2.86  , -2.79  , -2.72  , -2.65  , -2.62  ,
                   -2.59  , -2.56  , -2.53  , -2.5   , -2.4644, -2.4288, -2.3932,
                   -2.3576, -2.322 , -2.292 , -2.262 , -2.232 , -2.202 , -2.172 ,
                   -2.1462, -2.1204, -2.0946, -2.0688, -2.043 , -2.0204, -1.9978,
                   -1.9752, -1.9526, -1.93  , -1.91  , -1.89  , -1.87  , -1.85  ,
                   -1.83  , -1.812 , -1.794 , -1.776 , -1.758 , -1.74  , -1.7226,
                   -1.7052, -1.6878, -1.6704, -1.653 , -1.6356, -1.6182, -1.6008,
                   -1.5834, -1.566 , -1.5486, -1.5312, -1.5138, -1.4964, -1.479 }; // 56
// EE+(phi,eta) projection eta edges
static const double eta_bins_EEp[5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB)+1] =
                   {1.479 ,  1.4964,  1.5138,  1.5312,  1.5486,  1.566 ,  1.5834,
                    1.6008,  1.6182,  1.6356,  1.653 ,  1.6704,  1.6878,  1.7052,
                    1.7226,  1.74  ,  1.758 ,  1.776 ,  1.794 ,  1.812 ,  1.83  ,
                    1.85  ,  1.87  ,  1.89  ,  1.91  ,  1.93  ,  1.9526,  1.9752,
                    1.9978,  2.0204,  2.043 ,  2.0688,  2.0946,  2.1204,  2.1462,
                    2.172 ,  2.202 ,  2.232 ,  2.262 ,  2.292 ,  2.322 ,  2.3576,
                    2.3932,  2.4288,  2.4644,  2.5   ,  2.53  ,  2.56  ,  2.59  ,
                    2.62  ,  2.65  ,  2.72  ,  2.79  ,  2.86  ,  2.93  ,  3.    }; // 56
